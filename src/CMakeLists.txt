

###############################################################################
# FIND EXTERNAL LIBRARIES REQUIRED BY PERCOLATOR
###############################################################################

# INCLUDE HEADER FILES IN SRC DIRECTORY
include_directories(${CMAKE_CURRENT_BINARY_DIR})

MESSAGE( STATUS "Checking packages:" )

if(XML_SUPPORT)
  find_package(XercesC REQUIRED) # XSD DEPENDS ON ${XERCESC_FOUND}
  if(XERCESC_FOUND)
    message(STATUS "Found Xerces-C: ${XERCESC_LIBRARIES}")
    message(STATUS "              : ${XERCESC_INCLUDE_DIR}")
    message(STATUS "       Version: ${XERCESC_VERSION}")
  else(XERCESC_FOUND)
    message(FATAL_ERROR "The package Apache Xerces-C++ has not been found")
  endif(XERCESC_FOUND)
  include_directories(${XERCESC_INCLUDE_DIR})

  find_package(Threads REQUIRED) # XERCESC depends on pthread

  if(APPLE OR MINGW)
    find_package(CURL REQUIRED) # XERCESC depends on curl
  endif()

  message(STATUS "Locating CodeSynthesis library with XSDDIR set to $ENV{XSDDIR}")
  find_package(Xsd REQUIRED) # CODESYNTHESIS
  if(XSD_FOUND)
    message(STATUS "Package CodeSynthesis XSD found: ${XSD_INCLUDE_DIR}, ${XSD_EXECUTABLE}")
  else(XSD_FOUND)
    message(FATAL_ERROR "The package CodeSynthesis XSD has not been found")
  endif(XSD_FOUND)
  include_directories(${XSD_INCLUDE_DIR})

  if(UNIX)
      set(CMAKE_PREFIX_PATH ${RPCDIR})
      find_package(RPC QUIET)
      set(HAVE_RPC ${RPC_FOUND})
      if(RPC_FOUND)
          message(STATUS "RPC found, will include headers from: ${RPC_INCLUDE_DIR}")
          message(STATUS "and libraries from                  : ${RPC_LIBRARIES}")
          set(LIBRARIES ${LIBRARIES} ${RPC_LIBRARIES})
          include_directories(${RPC_INCLUDE_DIR})
      else(RPC_FOUND)
          message(FATAL_ERROR "RPC support is mandatory.\n"
          "Note that SunRPC has been removed in glibc-2.26 and later, "
          "while being optional in earlier versions. Consider using the "
          "recommended and more modern libtirpc instead.\n"
          "Use -DRPCDIR=DIR to specify the rpc directory tree.\n")
      endif(RPC_FOUND)
  endif()
endif(XML_SUPPORT)

if(APPLE AND OPENMP_FOUND)
  include_directories(${OpenMP_CXX_INCLUDE_DIRS})
endif(APPLE AND OPENMP_FOUND)

# Boost header libraries only for async io header library (asio.hpp)
find_package(Boost ${BOOST_MIN_VERSION} COMPONENTS filesystem system REQUIRED)
add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED)
add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)
if(WIN32)
  add_definitions(-DBOOST_ALL_NO_LIB) # disable autolinking in boost
endif(WIN32)
include_directories(${Boost_INCLUDE_DIRS})

#########################################
# COMPILE BLAS
#########################################
add_subdirectory(dblas)
set (DBLAS_LIBRARIES ${DBLAS_LIBRARIES} dblas)

###############################################################################
# PULL Eigen
###############################################################################

include(FetchContent)

# Avoid Eigen's conflicting 'blas' target
set(EIGEN_BUILD_BLAS OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "" FORCE)

# Avoid adding Eigen's blas directory (explicitly ignore it)
set(EIGEN_DISABLE_BLAS ON CACHE BOOL "Disable Eigen's BLAS build" FORCE)

FetchContent_Declare(
  eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3.4.0
)

FetchContent_MakeAvailable(eigen)
include_directories(${eigen_SOURCE_DIR})

### Patch for avoiding _finite

set(EIGEN_MATHFUNCTIONS_FILE "${eigen_SOURCE_DIR}/Eigen/src/Core/MathFunctions.h")

# Do the patch if not already patched
file(READ "${EIGEN_MATHFUNCTIONS_FILE}" MATHFUNCS_CONTENT)

string(REGEX REPLACE
  "#elif EIGEN_USE_STD_FPCLASSIFY[ \r\n]+[^\r\n]*using std::isfinite;[ \r\n]+[^\r\n]*return isfinite[^\r\n]*;"
  "#elif EIGEN_USE_STD_FPCLASSIFY\n    return !std::isnan(x) && !std::isinf(x); // PATCHED to avoid std::isfinite on MSVC"
  MATHFUNCS_CONTENT_PATCHED
  "${MATHFUNCS_CONTENT}"
)

file(WRITE "${EIGEN_MATHFUNCTIONS_FILE}" "${MATHFUNCS_CONTENT_PATCHED}")
### End patch

###############################################################################
# RUN CODESYNTHESIS
###############################################################################

if(XML_SUPPORT)
  add_definitions(-DXSD_CXX11)
  if(MINGW OR WIN32)
    set(xdr_flags)
  else()
    set(xdr_flags --generate-insertion XDR --generate-extraction XDR --std c++11)
  endif()

  # Generate percolator_in
  set(xsdfiles_in ${CMAKE_CURRENT_BINARY_DIR}/percolator_in.cxx ${CMAKE_CURRENT_BINARY_DIR}/percolator_in.hxx)
  set(percolator_in xml/percolator_in.xsd)
  add_custom_command(
    OUTPUT ${xsdfiles_in}
    COMMAND ${XSD_EXECUTABLE} cxx-tree --generate-serialization ${xdr_flags} --root-element-all
            --namespace-map ${percolator-in-namespace}=percolatorInNs
            ${CMAKE_CURRENT_SOURCE_DIR}/${percolator_in}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${percolator_in}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating percolator_in.* from XSD"
  )

  # Generate percolator_out
  set(xsdfiles_out ${CMAKE_CURRENT_BINARY_DIR}/percolator_out.cxx ${CMAKE_CURRENT_BINARY_DIR}/percolator_out.hxx)
  set(percolator_out xml/percolator_out.xsd)
  add_custom_command(
    OUTPUT ${xsdfiles_out}
    COMMAND ${XSD_EXECUTABLE} cxx-tree --generate-serialization ${xdr_flags} --root-element-all
            --namespace-map ${percolator-out-namespace}=percolatorOutNs
            ${CMAKE_CURRENT_SOURCE_DIR}/${percolator_out}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${percolator_out}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating percolator_out.* from XSD"
  )

  # 👇 This is the critical part!
  add_custom_target(generate_xsd ALL DEPENDS ${xsdfiles_in} ${xsdfiles_out})
endif()

#############################################################################
# COMPILE PERCOLATOR
###############################################################################

if(XML_SUPPORT)
  add_library(perclibrary STATIC ${xsdfiles_in} ${xsdfiles_out} parser.cxx serializer.cxx BaseSpline.cpp MassHandler.cpp
                  PSMDescription.cpp ResultHolder.cpp IsotonicPEP.cpp
								  XMLInterface.cpp SetHandler.cpp StdvNormalizer.cpp svm.cpp Caller.cpp CompositionSorter.cpp CrossValidation.cpp Enzyme.cpp Globals.cpp NoNormalizer.cpp Normalizer.cpp
								  SanityCheck.cpp UniNormalizer.cpp DataSet.cpp TabFileValidator.cpp FeatureNames.cpp LogisticRegression.cpp Numerical.cpp Option.cpp PosteriorEstimator.cpp
								  ProteinProbEstimator.cpp ProteinFDRestimator.cpp ScoreHolder.cpp Scores.cpp PseudoRandom.cpp Set.cpp SqtSanityCheck.cpp ssl.cpp PackedVector.cpp
								  PackedMatrix.cpp Logger.cpp MyException.cpp ProteinScoreHolder.cpp PickedProteinInterface.cpp FeatureMemoryPool.cpp Timer.cpp TmpDir.cpp ValidateTabFile.cpp Vector.cpp)
  add_dependencies(perclibrary generate_xsd)
else(XML_SUPPORT)
  add_library(perclibrary STATIC BaseSpline.cpp MassHandler.cpp ResultHolder.cpp PSMDescription.cpp IsotonicPEP.cpp
								  XMLInterface.cpp SetHandler.cpp StdvNormalizer.cpp svm.cpp Caller.cpp CompositionSorter.cpp CrossValidation.cpp Enzyme.cpp Globals.cpp NoNormalizer.cpp Normalizer.cpp
								  SanityCheck.cpp UniNormalizer.cpp DataSet.cpp TabFileValidator.cpp FeatureNames.cpp LogisticRegression.cpp Numerical.cpp Option.cpp PosteriorEstimator.cpp
								  ProteinProbEstimator.cpp ProteinFDRestimator.cpp ScoreHolder.cpp Scores.cpp PseudoRandom.cpp Set.cpp SqtSanityCheck.cpp ssl.cpp PackedVector.cpp
								  PackedMatrix.cpp Logger.cpp MyException.cpp ProteinScoreHolder.cpp PickedProteinInterface.cpp FeatureMemoryPool.cpp Timer.cpp TmpDir.cpp ValidateTabFile.cpp Vector.cpp)
endif(XML_SUPPORT)


###############################################################################
# COMPILE INTERNAL LIBRARIES
###############################################################################

#COMPILE PICKED PROTEIN
add_subdirectory(picked_protein)

#INCLUDE PERCOLATOR HEADERS
include_directories(.)

#INCLUDE PICKED PROTEIN HEADERS FOR PERCLIBRARY
include_directories(picked_protein)

################################################################################

# COMPILE PERCOLATOR
if(CRUX)
  ADD_DEFINITIONS(-DCRUX)
  include_directories(${CRUX} ${CRUX}/src ${EXT_BINARY_DIR})
ELSE(CRUX)
  add_executable(percolator main.cpp)

  if(APPLE)
    set_property(TARGET percolator PROPERTY LINK_SEARCH_START_STATIC FALSE)
    set_property(TARGET percolator PROPERTY LINK_SEARCH_END_STATIC TRUE)
  elseif(MINGW)
    set_target_properties(percolator PROPERTIES LINK_FLAGS "-Wl,-Bdynamic -static-libgcc")
  endif()

  if(PROFILING AND UNIX)
    set_target_properties(percolator PROPERTIES COMPILE_FLAGS -pg LINK_FLAGS -pg)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 ${PROFILING_COMMON_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 ${PROFILING_COMMON_FLAGS}")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 ${PROFILING_COMMON_FLAGS}")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 ${PROFILING_COMMON_FLAGS}")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -DNDEBUG ${PROFILING_COMMON_FLAGS}")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -g -DNDEBUG ${PROFILING_COMMON_FLAGS}")
  endif()

  # SET LIBRARIES FOR PERCOLATOR
  set(COMMON_LIBRARIES perclibrary ${DBLAS_LIBRARIES} picked_protein ${XERCESC_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${CURL_LIBRARIES} ${RPC_LIBRARIES})
  if(NOT MSVC)
    if(APPLE)
      set(COMMON_LIBRARIES ${COMMON_LIBRARIES} ${OpenMP_CXX_LIBRARIES})
    else(APPLE)
      set(COMMON_LIBRARIES ${COMMON_LIBRARIES} stdc++)
    endif(APPLE)
  endif(NOT MSVC)
  if(COVERAGE)
    target_link_libraries(percolator -fprofile-arcs ${COMMON_LIBRARIES} ${Boost_LIBRARIES})
  else(COVERAGE)
    target_link_libraries(percolator ${COMMON_LIBRARIES} ${Boost_LIBRARIES})
  endif(COVERAGE)

  # INSTALL PERCOLATOR
  install(TARGETS percolator EXPORT PERCOLATOR DESTINATION bin) # Important to use relative path here (used by CPack)!
ENDIF(CRUX)

###############################################################################
# COMPILE QVALITY
###############################################################################

add_subdirectory(qvality)

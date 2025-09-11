

###############################################################################
# FIND EXTERNAL LIBRARIES REQUIRED BY PERCOLATOR
###############################################################################

# INCLUDE HEADER FILES IN SRC DIRECTORY
include_directories(${CMAKE_CURRENT_BINARY_DIR})

MESSAGE( STATUS "Checking packages:" )

if(APPLE AND OPENMP_FOUND)
  include_directories(${OpenMP_CXX_INCLUDE_DIRS})
endif(APPLE AND OPENMP_FOUND)


if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 OLD)
endif()

# Link Boost static and single threaded
#if(APPLE)
#  set(Boost_USE_STATIC_LIBS OFF)
#  set(Boost_USE_MULTITHREADED ON)
#else(APPLE)
  set(Boost_USE_STATIC_LIBS   ON)
  set(Boost_USE_MULTITHREADED OFF)
#endif()
set(Boost_NO_BOOST_CMAKE ON CACHE BOOL "" FORCE)

# Boost header libraries only for async io header library (asio.hpp)
find_package(Boost ${BOOST_MIN_VERSION} COMPONENTS filesystem system REQUIRED)
add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED)
add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)
if(WIN32)
  add_definitions(-DBOOST_ALL_NO_LIB) # disable autolinking in boost
endif(WIN32)
if(MSVC)
  add_definitions (/MP /D _CRT_SECURE_NO_WARNINGS /D _CRT_NONSTDC_NO_DEPRECATE)
  # CMake automatically configures VS to use multi threading,
  # if boost would use single threading, the auto-linking process of boost would fail.
  set(Boost_USE_MULTITHREADED ON)
endif(MSVC)
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

#############################################################################
# COMPILE PERCOLATOR
###############################################################################

add_library(perclibrary STATIC BaseSpline.cpp MassHandler.cpp ResultHolder.cpp PSMDescription.cpp IsotonicPEP.cpp
								  SetHandler.cpp StdvNormalizer.cpp svm.cpp Caller.cpp CompositionSorter.cpp CrossValidation.cpp Enzyme.cpp Globals.cpp NoNormalizer.cpp Normalizer.cpp
								  SanityCheck.cpp UniNormalizer.cpp DataSet.cpp TabFileValidator.cpp FeatureNames.cpp LogisticRegression.cpp Numerical.cpp Option.cpp PosteriorEstimator.cpp
								  ProteinProbEstimator.cpp ProteinFDRestimator.cpp ScoreHolder.cpp Scores.cpp PseudoRandom.cpp Set.cpp SqtSanityCheck.cpp ssl.cpp PackedVector.cpp
								  PackedMatrix.cpp Logger.cpp MyException.cpp ProteinScoreHolder.cpp PickedProteinInterface.cpp FeatureMemoryPool.cpp Timer.cpp TmpDir.cpp ValidateTabFile.cpp Vector.cpp)


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

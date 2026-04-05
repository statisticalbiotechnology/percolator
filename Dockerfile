FROM docker.io/library/ubuntu:24.04 AS builder
ARG percolator_cmake_args="-DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DGOOGLE_TEST=0"

RUN apt-get update && apt-get install -y 

# Remove the system gtest libraries (if installed)
RUN apt-get remove --purge -y libgtest-dev && apt-get autoremove -y

RUN apt-get install -y build-essential g++ make cmake gawk git ca-certificates

RUN apt-get install -y \
    libboost-filesystem-dev \
    libboost-system-dev \
    libboost-thread-dev 


COPY / /percolator

RUN mkdir -p /release /build
RUN mkdir -p /build/percolator;

WORKDIR /build/percolator

RUN cmake ${percolator_cmake_args} /percolator
RUN make
RUN make install

FROM docker.io/library/ubuntu:24.04 AS runtime

RUN apt-get update && apt-get install -y \
    libboost-system-dev \
    libboost-filesystem-dev \
    libboost-thread-dev

COPY --from=builder /usr/bin/percolator /usr/bin/percolator
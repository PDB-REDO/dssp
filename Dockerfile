FROM ubuntu:20.04

ENV TZ="America/New_York"
RUN apt-get update && \
    apt-get install -yq tzdata && \
    ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime && \
    dpkg-reconfigure -f noninteractive tzdata && \
    apt install -y build-essential cmake zlib1g-dev git libeigen3-dev

WORKDIR /build

# Build and install libcifpp
# https://github.com/PDB-REDO/libcifpp
RUN cd /build && \
    git clone https://github.com/PDB-REDO/libcifpp.git --recurse-submodules && \
    cd libcifpp && \
    cmake -S . -B build -DCMAKE_INSTALL_PREFIX=/root/.local -DCMAKE_BUILD_TYPE=Release && \
    cmake --build build && \
    cmake --install build && \
    echo "libcifpp installed"

# Build and install libmcfp
# https://github.com/mhekkel/libmcfp
RUN cd /build && \
    git clone https://github.com/mhekkel/libmcfp.git && \
    cd libmcfp && \
    mkdir build && \
    cd build && \
    cmake .. && \
    cmake --build . && \
    cmake --install . && \
    echo "libmcfp installed"

# Build and install dssp
COPY . /src
RUN cd /src && \
    mkdir build && \
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/root/.local/lib/cmake/cifpp/ && \
    cmake --build build && \
    cmake --install build && \
    echo "dssp installed" && \
    rm -rf /src /build

WORKDIR /data
ENTRYPOINT ["mkdssp"]
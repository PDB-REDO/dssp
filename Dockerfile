FROM ubuntu:22.04

ENV TZ="Europe/Amsterdam"
RUN apt-get update && \
    apt-get install -yq tzdata && \
    ln -fs /usr/share/zoneinfo/Europe/Amsterdam /etc/localtime && \
    dpkg-reconfigure -f noninteractive tzdata && \
    apt install -y build-essential cmake zlib1g-dev git libeigen3-dev

WORKDIR /build

# Build and install libcifpp
# https://github.com/PDB-REDO/libcifpp
RUN cd /build && \
    git clone https://github.com/PDB-REDO/libcifpp.git && \
    cd libcifpp && \
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
		-DBUILD_TESTING=OFF -DCIFPP_DOWNLOAD_CCD=OFF && \
    cmake --build build -j $(nproc) && \
    cmake --install build && \
    echo "libcifpp installed"

# Build and install libmcfp
# https://github.com/mhekkel/libmcfp
RUN cd /build && \
    git clone https://github.com/mhekkel/libmcfp.git && \
    cd libmcfp && \
    cmake -S . -B build -DBUILD_TESTING=OFF && \
    cmake --build build -j $(nproc) && \
    cmake --install build && \
    echo "libmcfp installed"

# Build and install dssp
COPY . /src
RUN cd /src && \
	rm -rf build && \
    mkdir build && \
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF && \
    cmake --build build -j $(nproc) && \
    cmake --install build && \
    echo "dssp installed" && \
    rm -rf /src /build

WORKDIR /data
ENTRYPOINT ["mkdssp"]
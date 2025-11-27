FROM ubuntu:24.04

ENV TZ="Europe/Amsterdam"
RUN apt-get update && \
    apt-get install -yq tzdata && \
    ln -fs /usr/share/zoneinfo/Europe/Amsterdam /etc/localtime && \
    dpkg-reconfigure -f noninteractive tzdata && \
    apt install -y build-essential cmake zlib1g-dev git libeigen3-dev libpcre2-dev

WORKDIR /build

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

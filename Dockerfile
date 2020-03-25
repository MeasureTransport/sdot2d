FROM debian:latest

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        binutils \
        ca-certificates \
        cmake \
        g++ \
        gcc \
        git \
        wget \
        build-essential \
        libgmp-dev \
        libcgal-dev \
    && rm -rf /var/lib/apt/lists/*

RUN useradd -rm -d /home/cgal -s /bin/bash cgal
USER cgal
WORKDIR /home/cgal

FROM jupyter/scipy-notebook

USER root

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

USER $NB_UID
WORKDIR /home/sdot

FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && apt-get install -y \
  build-essential \
  curl \
  wget \
  git \
  flex \
  bison \
  texinfo \
  libgmp-dev \
  libmpfr-dev \
  libmpc-dev \
  libisl-dev \
  zlib1g-dev \
  python3 \
  cmake \
  ninja-build \
  unzip \
  patch \
  pkg-config \
  xz-utils \
  libtclap-dev \
  libopenvdb-dev \
  libtbb-dev \
  libcppunit-dev \
  libeigen3-dev \
  liblzma-dev \
  zlib1g-dev \
  libbz2-dev \
  libssl-dev

RUN apt-get install curl wget

# -------- Build & install latest Boost --------
ARG BOOST_PREFIX="/usr/local"
ARG BOOST_WITH="--with-system --with-filesystem --with-thread --with-chrono --with-date_time --with-regex --with-program_options --with-iostreams"

WORKDIR /tmp
RUN wget https://archives.boost.io/release/1.89.0/source/boost_1_89_0.tar.bz2
RUN tar -xvjf boost_1_89_0.tar.bz2
WORKDIR /tmp/boost_1_89_0
RUN ./bootstrap.sh --prefix="${BOOST_PREFIX}"
RUN ./b2 -j"$(nproc)" ${BOOST_WITH} variant=release runtime-link=shared install
RUN echo "${BOOST_PREFIX}/lib" > /etc/ld.so.conf.d/boost.conf
RUN ldconfig

# -------- Cleanup --------
RUN apt-get clean

CMD ["g++", "--version"]
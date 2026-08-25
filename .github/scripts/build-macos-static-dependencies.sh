#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <architecture> <install-prefix>" >&2
    exit 2
fi

architecture="$1"
install_prefix="$2"
mkdir -p "$install_prefix"
install_prefix="$(cd "$install_prefix" && pwd)"

deployment_target="12.0"
sdk_path="$(xcrun --sdk macosx --show-sdk-path)"
build_jobs="$(sysctl -n hw.logicalcpu)"
source_dir="$(mktemp -d)"
trap 'rm -rf "$source_dir"' EXIT

export MACOSX_DEPLOYMENT_TARGET="$deployment_target"
export CFLAGS="-arch $architecture -isysroot $sdk_path -mmacosx-version-min=$deployment_target"
export CXXFLAGS="$CFLAGS"
export LDFLAGS="$CFLAGS"
export CPATH="$install_prefix/include${CPATH:+:$CPATH}"
export LIBRARY_PATH="$install_prefix/lib${LIBRARY_PATH:+:$LIBRARY_PATH}"

download() {
    local filename="$1"
    local checksum="$2"
    local url="$3"

    curl --fail --location --retry 3 --output "$filename" "$url"
    echo "$checksum  $filename" | shasum -a 256 --check
}

cd "$source_dir"

download xz-5.8.3.tar.gz \
    3d3a1b973af218114f4f889bbaa2f4c037deaae0c8e815eec381c3d546b974a0 \
    https://github.com/tukaani-project/xz/releases/download/v5.8.3/xz-5.8.3.tar.gz
download zstd-1.5.7.tar.gz \
    37d7284556b20954e56e1ca85b80226768902e2edabd3b649e9e72c0c9012ee3 \
    https://github.com/facebook/zstd/archive/refs/tags/v1.5.7.tar.gz
download lz4-1.10.0.tar.gz \
    537512904744b35e232912055ccf8ec66d768639ff3abe5788d90d792ec5f48b \
    https://github.com/lz4/lz4/archive/refs/tags/v1.10.0.tar.gz
download snappy-1.2.2.tar.gz \
    90f74bc1fbf78a6c56b3c4a082a05103b3a56bb17bca1a27e052ea11723292dc \
    https://github.com/google/snappy/archive/refs/tags/1.2.2.tar.gz
download blosc-1.21.6.tar.gz \
    9fcd60301aae28f97f1301b735f966cc19e7c49b6b4321b839b4579a0c156f38 \
    https://github.com/Blosc/c-blosc/archive/refs/tags/v1.21.6.tar.gz
download boost-1.90.0.tar.xz \
    9e6bee9ab529fb2b0733049692d57d10a72202af085e553539a05b4204211a6f \
    https://github.com/boostorg/boost/releases/download/boost-1.90.0/boost-1.90.0-b2-nodocs.tar.xz

tar -xzf xz-5.8.3.tar.gz
tar -xzf zstd-1.5.7.tar.gz
tar -xzf lz4-1.10.0.tar.gz
tar -xzf snappy-1.2.2.tar.gz
tar -xzf blosc-1.21.6.tar.gz
tar -xJf boost-1.90.0.tar.xz

(
    cd xz-5.8.3
    ./configure \
        --prefix="$install_prefix" \
        --disable-shared \
        --enable-static \
        --disable-doc
    make -j"$build_jobs"
    make install
)

cmake -S zstd-1.5.7/build/cmake -B zstd-build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$install_prefix" \
    -DCMAKE_INSTALL_LIBDIR=lib \
    -DCMAKE_OSX_ARCHITECTURES="$architecture" \
    -DCMAKE_OSX_DEPLOYMENT_TARGET="$deployment_target" \
    -DZSTD_BUILD_PROGRAMS=OFF \
    -DZSTD_BUILD_SHARED=OFF \
    -DZSTD_BUILD_STATIC=ON \
    -DZSTD_BUILD_TESTS=OFF
cmake --build zstd-build --parallel "$build_jobs"
cmake --install zstd-build

cmake -S lz4-1.10.0/build/cmake -B lz4-build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$install_prefix" \
    -DCMAKE_INSTALL_LIBDIR=lib \
    -DCMAKE_OSX_ARCHITECTURES="$architecture" \
    -DCMAKE_OSX_DEPLOYMENT_TARGET="$deployment_target" \
    -DBUILD_SHARED_LIBS=OFF \
    -DLZ4_BUILD_CLI=OFF \
    -DLZ4_BUILD_LEGACY_LZ4C=OFF
cmake --build lz4-build --parallel "$build_jobs"
cmake --install lz4-build

cmake -S snappy-1.2.2 -B snappy-build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$install_prefix" \
    -DCMAKE_INSTALL_LIBDIR=lib \
    -DCMAKE_OSX_ARCHITECTURES="$architecture" \
    -DCMAKE_OSX_DEPLOYMENT_TARGET="$deployment_target" \
    -DBUILD_SHARED_LIBS=OFF \
    -DSNAPPY_BUILD_BENCHMARKS=OFF \
    -DSNAPPY_BUILD_TESTS=OFF
cmake --build snappy-build --parallel "$build_jobs"
cmake --install snappy-build

cmake -S c-blosc-1.21.6 -B blosc-build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
    -DCMAKE_INSTALL_PREFIX="$install_prefix" \
    -DCMAKE_INSTALL_LIBDIR=lib \
    -DCMAKE_OSX_ARCHITECTURES="$architecture" \
    -DCMAKE_OSX_DEPLOYMENT_TARGET="$deployment_target" \
    -DCMAKE_PREFIX_PATH="$install_prefix" \
    -DBUILD_STATIC=ON \
    -DBUILD_SHARED=OFF \
    -DBUILD_TESTS=OFF \
    -DBUILD_BENCHMARKS=OFF \
    -DBUILD_FUZZERS=OFF \
    -DDEACTIVATE_AVX2=ON \
    -DDEACTIVATE_SNAPPY=ON \
    -DPREFER_EXTERNAL_LZ4=ON \
    -DPREFER_EXTERNAL_ZLIB=ON \
    -DPREFER_EXTERNAL_ZSTD=ON
cmake --build blosc-build --parallel "$build_jobs"
cmake --install blosc-build

(
    cd boost-1.90.0
    ./bootstrap.sh \
        --prefix="$install_prefix" \
        --with-libraries=filesystem,iostreams,regex \
        --without-icu
    ./b2 -j"$build_jobs" install \
        --prefix="$install_prefix" \
        --libdir="$install_prefix/lib" \
        --layout=system \
        variant=release \
        link=static \
        runtime-link=shared \
        threading=multi \
        "cxxflags=$CXXFLAGS -I$install_prefix/include" \
        "linkflags=$LDFLAGS -L$install_prefix/lib"
)

required_archives=(
    libblosc.a
    libboost_filesystem.a
    libboost_iostreams.a
    libboost_regex.a
    liblz4.a
    liblzma.a
    libsnappy.a
    libzstd.a
)
for archive in "${required_archives[@]}"; do
    if [[ ! -f "$install_prefix/lib/$archive" ]]; then
        echo "Static dependency was not produced: $install_prefix/lib/$archive" >&2
        exit 1
    fi
    if ! lipo -verify_arch "$architecture" "$install_prefix/lib/$archive"; then
        echo "$archive does not contain the required $architecture architecture." >&2
        exit 1
    fi
done

cat > "$install_prefix/VERSIONS.txt" <<'EOF'
xz 5.8.3
zstd 1.5.7
lz4 1.10.0
snappy 1.2.2
c-blosc 1.21.6
boost 1.90.0
minimum macOS 12.0
EOF

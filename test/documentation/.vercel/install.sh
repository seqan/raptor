#!/usr/bin/env bash
set -exo pipefail

DOXYGEN_VERSION=1.9.5
WORK_DIR=`pwd`
CACHE_DIR="${WORK_DIR}/node_modules"

# Install dependencies.
amazon-linux-extras install epel &>/dev/null
yum --assumeyes --quiet install wget cmake3 flex bison xz graphviz xapian-core-devel &>/dev/null

# Download doxygen.
mkdir -p ${CACHE_DIR}/doxygen-download
wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --quiet --directory-prefix=${CACHE_DIR}/doxygen-download/ https://sourceforge.net/projects/doxygen/files/rel-${DOXYGEN_VERSION}/doxygen-${DOXYGEN_VERSION}.src.tar.gz
tar -C ${CACHE_DIR}/doxygen-download -zxf ${CACHE_DIR}/doxygen-download/doxygen-${DOXYGEN_VERSION}.src.tar.gz

# Configure doxygen.
cd ${CACHE_DIR}/doxygen-download/doxygen-${DOXYGEN_VERSION}
mkdir -p build && cd build
cmake_command() {
    cmake3 -Dbuild_search=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-D_GLIBCXX_USE_CXX11_ABI=0 -w" -G "Unix Makefiles" .. 1>/dev/null
}
if ! cmake_command; then
    rm -rf *
    cmake_command
fi

# Build doxygen.
make -j 4 1>/dev/null

# Install doxygen
make install DESTDIR=${CACHE_DIR}/doxygen-${DOXYGEN_VERSION} 1>/dev/null

# Symlink doxygen.
ln -s ${CACHE_DIR}/doxygen-${DOXYGEN_VERSION}/usr/local/bin ${WORK_DIR}/doxygen-bin
cp /usr/lib64/libxapian.so ${WORK_DIR}/doxygen-bin

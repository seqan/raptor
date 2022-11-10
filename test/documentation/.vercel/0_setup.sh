#!/usr/bin/env bash
set -euxo pipefail

DOXYGEN_VERSION=1.9.5
WORK_DIR=`pwd`
CACHE_DIR="${WORK_DIR}/node_modules" # The node_modules directory is always cached.

### Install dependencies.
# Enables EPEL 7 which hosts xapian-core-devel. The server runs on Amazon Linux 2.
amazon-linux-extras install epel &>/dev/null
yum --assumeyes --quiet install wget cmake3 flex bison xz graphviz xapian-core-devel &>/dev/null

### Download doxygen.
mkdir -p ${CACHE_DIR}/doxygen-download
wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --quiet --directory-prefix=${CACHE_DIR}/doxygen-download/ https://sourceforge.net/projects/doxygen/files/rel-${DOXYGEN_VERSION}/doxygen-${DOXYGEN_VERSION}.src.tar.gz
tar -C ${CACHE_DIR}/doxygen-download -zxf ${CACHE_DIR}/doxygen-download/doxygen-${DOXYGEN_VERSION}.src.tar.gz

### Configure doxygen.
cd ${CACHE_DIR}/doxygen-download/doxygen-${DOXYGEN_VERSION}
mkdir -p build && cd build

# xapian-core-devel from epel (EPEL 7) is built with pre-cxx11-abi, so we need to do so for doxygen, too.
# -Dbuild_search=ON enables building doxyindexer and doxysearch.cgi
cmake_command() {
    cmake3 .. \
        -Dbuild_search=ON \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_FLAGS="-D_GLIBCXX_USE_CXX11_ABI=0 -w" \
        -G "Unix Makefiles" 1>/dev/null
}

# If CMake fails, we clear the directory and try again. Might happen due to path changes on vercel's side.
if ! cmake_command; then
    rm -rf *
    cmake_command
fi

### Build doxygen.
make -j4 1>/dev/null

### Install doxygen.
make install DESTDIR=${CACHE_DIR}/doxygen-${DOXYGEN_VERSION} 1>/dev/null

### Symlink.
# Just allows using an easier path.
ln -s ${CACHE_DIR}/doxygen-${DOXYGEN_VERSION}/usr/local/bin ${WORK_DIR}/doxygen-bin
# We will need to ship the dynamically linked libxapian. Copied for easy access in api/doxysearch.sh
cp /usr/lib64/libxapian.so ${WORK_DIR}/doxygen-bin

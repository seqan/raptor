#!/usr/bin/env bash
set -exo pipefail

SOURCE_DIR=`pwd`
CACHE_DIR="${SOURCE_DIR}/node_modules"
BUILD_DIR="${CACHE_DIR}/build"
EXPORT_DIR="${SOURCE_DIR}/export"
PATH="$PATH:${SOURCE_DIR}/doxygen-bin"

# Print versions.
${SOURCE_DIR}/doxygen-bin/doxygen --version
cmake3 --version
doxygen --version

# Configure documentation build.
mkdir -p "${BUILD_DIR}" && cd "${BUILD_DIR}"
cmake3 -DRAPTOR_VERCEL_PREVIEW_DOC=ON -DCMAKE_INSTALL_PREFIX="" -DCMAKE_INSTALL_DOCDIR="." "${SOURCE_DIR}/test/documentation" 1>/dev/null

# Build documentation.
cmake3 --build . --target download-cppreference-doxygen-web-tag 1>/dev/null
make doc &> /dev/null
# ctest3 -j 2 --output-on-failure .

# Install documentation.
mkdir -p "${EXPORT_DIR}/"
DESTDIR="${EXPORT_DIR}/" cmake3 -DCOMPONENT=doc -P cmake_install.cmake 1>/dev/null

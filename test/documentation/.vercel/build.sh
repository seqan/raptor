#!/usr/bin/env bash
set -exo pipefail

WORK_DIR=`pwd`
CACHE_DIR="${WORK_DIR}/node_modules"
BUILD_DIR="${CACHE_DIR}/build"
EXPORT_DIR="${WORK_DIR}/export"
PATH="$PATH:${WORK_DIR}/doxygen-bin"

# Print versions.
cmake3 --version
doxygen --version

# Configure documentation build.
mkdir -p "${BUILD_DIR}" && cd "${BUILD_DIR}"
cmake3 -DRAPTOR_VERCEL_PREVIEW_DOC=ON -DRAPTOR_VERCEL_URL="${VERCEL_URL}" -DCMAKE_INSTALL_PREFIX="" -DCMAKE_INSTALL_DOCDIR="." "${WORK_DIR}/.." 1>/dev/null

# Build documentation.
cmake3 --build . --target download-cppreference-doxygen-web-tag 1>/dev/null
make doc &> /dev/null

# Install documentation.
mkdir -p "${EXPORT_DIR}/"
DESTDIR="${EXPORT_DIR}/" cmake3 -DCOMPONENT=doc -P cmake_install.cmake 1>/dev/null
cd ${WORK_DIR}/doxygen-bin
doxyindexer "${BUILD_DIR}/searchdata.xml"

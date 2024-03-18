#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

set -euxo pipefail

WORK_DIR=`pwd`
EXPORT_DIR="${WORK_DIR}/export"      # ${EXPORT_DIR}/html will be served by vercel.
PATH="${WORK_DIR}/doxygen-bin:$PATH" # Add binaries from 0_setup.sh to PATH.
CACHE_DIR="${WORK_DIR}/node_modules" # The node_modules directory is always cached.
BUILD_DIR="${CACHE_DIR}/build"       # We also want to cache the documentation build.

### Print versions.
cmake3 --version
doxygen --version

### Configure documentation build.
mkdir -p "${BUILD_DIR}" && cd "${BUILD_DIR}"
cmake3 "${WORK_DIR}/.." \
    -DRAPTOR_VERCEL_PREVIEW_DOC=ON \
    -DRAPTOR_VERCEL_URL="${VERCEL_URL}" \
    -DCMAKE_INSTALL_PREFIX="" \
    -DCMAKE_INSTALL_DOCDIR="."  1>/dev/null

### Build documentation.
cmake3 --build . --target download-cppreference-doxygen-web-tag 1>/dev/null
make -j4 doc &> /dev/null

### Install documentation.
mkdir -p "${EXPORT_DIR}/"
DESTDIR="${EXPORT_DIR}/" cmake3 -DCOMPONENT=doc -P cmake_install.cmake 1>/dev/null

### Run indexer.
# We want the resulting index to be in the binary dir for easy access in api/doxysearch.sh.
cd ${WORK_DIR}/doxygen-bin
# Will put a directory doxysearch.db in the current directory.
doxyindexer "${BUILD_DIR}/searchdata.xml"

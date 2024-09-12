# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

# CPM Package Lock
# This file should be committed to version control

# hibf
set (RAPTOR_HIBF_VERSION 485b22b3a770e35a6e8a143afae91024025c163a)
CPMDeclarePackage (hibf
                   NAME hibf
                   GIT_TAG ${RAPTOR_HIBF_VERSION}
                   GITHUB_REPOSITORY seqan/hibf
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "INSTALL_HIBF OFF"
)
# sharg
set (RAPTOR_SHARG_VERSION 7534545c5f5dca3876ca3d3b07c369131b185f92)
CPMDeclarePackage (sharg
                   NAME sharg
                   GIT_TAG ${RAPTOR_SHARG_VERSION}
                   GITHUB_REPOSITORY seqan/sharg-parser
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "INSTALL_SHARG OFF" "INSTALL_TDL OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING" "BUILD_TESTING OFF"
)
# seqan3
set (RAPTOR_SEQAN3_VERSION abcad130a70b86f2286ed0d80ef7bf3de3fc1e39)
CPMDeclarePackage (seqan3
                   NAME seqan3
                   GIT_TAG ${RAPTOR_SEQAN3_VERSION}
                   GITHUB_REPOSITORY seqan/seqan3
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "INSTALL_SEQAN3 OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)
# chopper
set (RAPTOR_CHOPPER_VERSION fe0e76307c77e8773dbdff2a4a56a0723dfd1ebe)
CPMDeclarePackage (chopper
                   NAME chopper
                   GIT_TAG ${RAPTOR_CHOPPER_VERSION}
                   GITHUB_REPOSITORY seqan/chopper
                   SYSTEM TRUE
                   OPTIONS "CHOPPER_INSTALL OFF" "CHOPPER_BUILD_DOC OFF" "CHOPPER_BUILD_TEST OFF"
                           "CMAKE_MESSAGE_LOG_LEVEL WARNING"
                   EXCLUDE_FROM_ALL TRUE
)
# benchmark
set (RAPTOR_BENCHMARK_VERSION 1.9.0)
CPMDeclarePackage (benchmark
                   NAME benchmark
                   VERSION ${RAPTOR_BENCHMARK_VERSION}
                   GITHUB_REPOSITORY google/benchmark
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "BENCHMARK_ENABLE_TESTING OFF" "BENCHMARK_ENABLE_WERROR OFF"
                           "CMAKE_MESSAGE_LOG_LEVEL WARNING" "CMAKE_CXX_STANDARD ${CMAKE_CXX_STANDARD}"
)
# googletest
set (RAPTOR_GOOGLETEST_VERSION 1.15.2)
CPMDeclarePackage (googletest
                   NAME googletest
                   VERSION ${RAPTOR_GOOGLETEST_VERSION}
                   GITHUB_REPOSITORY google/googletest
                   SYSTEM TRUE
                   OPTIONS "BUILD_GMOCK OFF" "INSTALL_GTEST OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING"
                           "CMAKE_CXX_STANDARD ${CMAKE_CXX_STANDARD}"
)
# use_ccache
set (USE_CCACHE_VERSION d2a54ef555b6fc2d496a4c9506dbeb7cf899ce37)
CPMDeclarePackage (use_ccache
                   NAME use_ccache
                   GIT_TAG ${USE_CCACHE_VERSION}
                   GITHUB_REPOSITORY seqan/cmake-scripts
                   SOURCE_SUBDIR ccache
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
)
# doxygen-awesome
set (RAPTOR_DOXYGEN_AWESOME_VERSION 2.3.3)
CPMDeclarePackage (doxygen_awesome
                   NAME doxygen_awesome
                   VERSION ${RAPTOR_DOXYGEN_AWESOME_VERSION}
                   GITHUB_REPOSITORY jothepro/doxygen-awesome-css
                   DOWNLOAD_ONLY TRUE
)

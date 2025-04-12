# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

# CPM Package Lock
# This file should be committed to version control

# cmake-format: off

# hibf
set (RAPTOR_HIBF_VERSION 5633fc044e93020984ad045c7001ede453b4b9e0 CACHE STRING "" FORCE)
CPMDeclarePackage (hibf
                   NAME hibf
                   GIT_TAG ${RAPTOR_HIBF_VERSION} # main
                   GITHUB_REPOSITORY seqan/hibf
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "INSTALL_HIBF OFF"
)
# sharg
set (RAPTOR_SHARG_VERSION dfff01056dda9271b158d34427f2d28fad9f7440 CACHE STRING "" FORCE)
CPMDeclarePackage (sharg
                   NAME sharg
                   GIT_TAG ${RAPTOR_SHARG_VERSION} # main
                   GITHUB_REPOSITORY seqan/sharg-parser
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "INSTALL_SHARG OFF" "INSTALL_TDL OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING" "BUILD_TESTING OFF"
)
# seqan3
set (RAPTOR_SEQAN3_VERSION 7f01aac8985b77b1363e6201316dc636b6d5f16d CACHE STRING "" FORCE)
CPMDeclarePackage (seqan3
                   NAME seqan3
                   GIT_TAG ${RAPTOR_SEQAN3_VERSION} # main
                   GITHUB_REPOSITORY seqan/seqan3
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "INSTALL_SEQAN3 OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)
# chopper
set (RAPTOR_CHOPPER_VERSION bef6e887af4e4f48262b1b192c68e863da26b688 CACHE STRING "" FORCE)
CPMDeclarePackage (chopper
                   NAME chopper
                   GIT_TAG ${RAPTOR_CHOPPER_VERSION} # main
                   GITHUB_REPOSITORY seqan/chopper
                   SYSTEM TRUE
                   OPTIONS "CHOPPER_INSTALL OFF" "CHOPPER_BUILD_DOC OFF" "CHOPPER_BUILD_TEST OFF"
                           "CMAKE_MESSAGE_LOG_LEVEL WARNING"
                   EXCLUDE_FROM_ALL TRUE
)
# benchmark
set (RAPTOR_BENCHMARK_VERSION 1.9.2 CACHE STRING "" FORCE)
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
set (RAPTOR_GOOGLETEST_VERSION 1.16.0 CACHE STRING "" FORCE)
CPMDeclarePackage (GTest
                   NAME GTest
                   VERSION ${RAPTOR_GOOGLETEST_VERSION}
                   GITHUB_REPOSITORY google/googletest
                   SYSTEM TRUE
                   OPTIONS "BUILD_GMOCK OFF" "INSTALL_GTEST OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING"
                           "CMAKE_CXX_STANDARD ${CMAKE_CXX_STANDARD}"
)
# use_ccache
set (USE_CCACHE_VERSION d2a54ef555b6fc2d496a4c9506dbeb7cf899ce37 CACHE STRING "" FORCE)
CPMDeclarePackage (use_ccache
                   NAME use_ccache
                   GIT_TAG ${USE_CCACHE_VERSION} # main
                   GITHUB_REPOSITORY seqan/cmake-scripts
                   SOURCE_SUBDIR ccache
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
)
# doxygen-awesome
set (RAPTOR_DOXYGEN_AWESOME_VERSION 2.3.4 CACHE STRING "" FORCE)
CPMDeclarePackage (doxygen_awesome
                   NAME doxygen_awesome
                   VERSION ${RAPTOR_DOXYGEN_AWESOME_VERSION}
                   GITHUB_REPOSITORY jothepro/doxygen-awesome-css
                   DOWNLOAD_ONLY TRUE
)
# ibf-fpga
set (RAPTOR_IBF_FPGA_VERSION 01533444a8f640cf043a369a855be082dd569dff CACHE STRING "" FORCE)
CPMDeclarePackage (ibf-fpga
                   NAME ibf-fpga
                   URL https://github.com/seqan/ibf-fpga/archive/${RAPTOR_IBF_FPGA_VERSION}.tar.gz # master
                   DOWNLOAD_ONLY YES
                   QUIET YES)

# cmake-format: on

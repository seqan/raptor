# SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

# CPM Package Lock
# This file should be committed to version control

# hibf
set (RAPTOR_HIBF_VERSION 2075a68d10f571399b5ff885f0b4cb2c94845601)
CPMDeclarePackage (hibf
                   NAME hibf
                   GIT_TAG ${RAPTOR_HIBF_VERSION}
                   GITHUB_REPOSITORY seqan/hibf
                   SYSTEM TRUE
                   OPTIONS "INSTALL_HIBF OFF"
)
# sharg
set (RAPTOR_SHARG_VERSION 3e96d37387125bd876ace3fdd11a20039ed0d21c)
CPMDeclarePackage (sharg
                   NAME sharg
                   GIT_TAG ${RAPTOR_SHARG_VERSION}
                   GITHUB_REPOSITORY seqan/sharg-parser
                   SYSTEM TRUE
                   OPTIONS "INSTALL_SHARG OFF" "INSTALL_TDL OFF"
)
# seqan3
set (RAPTOR_SEQAN3_VERSION 6268f1d11c58abf24c78c7f99a71817d13d5bddd)
CPMDeclarePackage (seqan3
                   NAME seqan3
                   GIT_TAG ${RAPTOR_SEQAN3_VERSION}
                   GITHUB_REPOSITORY seqan/seqan3
                   SYSTEM TRUE
                   OPTIONS "INSTALL_SEQAN3 OFF"
)
# chopper
set (RAPTOR_CHOPPER_VERSION bdcecb6d0f022778ae690da419654f26c7989d2c)
CPMDeclarePackage (chopper
                   NAME chopper
                   GIT_TAG ${RAPTOR_CHOPPER_VERSION}
                   GITHUB_REPOSITORY seqan/chopper
                   SYSTEM TRUE
                   OPTIONS "CHOPPER_INSTALL OFF" "CHOPPER_BUILD_DOC OFF" "CHOPPER_BUILD_TEST OFF"
                   EXCLUDE_FROM_ALL TRUE
)
# benchmark
set (RAPTOR_BENCHMARK_VERSION 1.8.2)
CPMDeclarePackage (benchmark
                   NAME benchmark
                   VERSION ${RAPTOR_BENCHMARK_VERSION}
                   GITHUB_REPOSITORY google/benchmark
                   SYSTEM TRUE
                   OPTIONS "BENCHMARK_ENABLE_TESTING OFF" "BENCHMARK_ENABLE_WERROR OFF"
)
# googletest
set (RAPTOR_GOOGLETEST_VERSION 1.14.0)
CPMDeclarePackage (googletest
                   NAME GTest
                   VERSION ${RAPTOR_GOOGLETEST_VERSION}
                   GITHUB_REPOSITORY google/googletest
                   SYSTEM TRUE
                   OPTIONS "BUILD_GMOCK OFF" "INSTALL_GTEST OFF"
)

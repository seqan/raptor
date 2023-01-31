# --------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
# --------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.10)

project (raptor_test_benchmark CXX)

set (RAPTOR_ENABLE_BENCHMARK TRUE)
include (../raptor-test.cmake)

raptor_require_benchmark ()
raptor_require_ccache ()
raptor_require_test ()

add_custom_target (cli_test)

add_cli_test (bin_influence_benchmark.cpp)

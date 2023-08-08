# --------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
# --------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.10)

project (raptor_test_performance LANGUAGES CXX C)

option (BUILD_RAPTOR_APP_FOR_TEST "" OFF)
include (../raptor-test.cmake)

raptor_require_benchmark ()
raptor_require_ccache ()

raptor_add_benchmark (bin_influence_benchmark.cpp)

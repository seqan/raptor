# --------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
# --------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.15)

add_cli_test (cli_build_hibf_chopper_test.cpp)
add_cli_test (cli_build_hibf_test.cpp)
add_cli_test (cli_build_ibf_compressed_test.cpp)
add_cli_test (cli_build_ibf_partitioned_test.cpp)
add_cli_test (cli_build_ibf_test.cpp)

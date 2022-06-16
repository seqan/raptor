# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.15)

project (raptor_utility_test)

get_filename_component(RAPTOR_UTIL_ROOT "${CMAKE_CURRENT_LIST_DIR}/../dummy" DIRECTORY)

add_subdirectory ("${RAPTOR_UTIL_ROOT}" "${CMAKE_CURRENT_BINARY_DIR}/util")
add_subdirectory ("${RAPTOR_UTIL_ROOT}/new" "${CMAKE_CURRENT_BINARY_DIR}/new")

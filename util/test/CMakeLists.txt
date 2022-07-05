# --------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
# --------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.15)

project (raptor_utility_test)

get_filename_component(RAPTOR_UTIL_ROOT "${CMAKE_CURRENT_LIST_DIR}/../dummy" DIRECTORY)

message (STATUS "Adding util...")
add_subdirectory ("${RAPTOR_UTIL_ROOT}" "${CMAKE_CURRENT_BINARY_DIR}/util")
message (STATUS "Adding new...")
add_subdirectory ("${RAPTOR_UTIL_ROOT}/new" "${CMAKE_CURRENT_BINARY_DIR}/new")
message (STATUS "Adding thresholding...")
add_subdirectory ("${RAPTOR_UTIL_ROOT}/thresholding" "${CMAKE_CURRENT_BINARY_DIR}/thresholding")
message (STATUS "Adding helper...")
add_subdirectory ("${RAPTOR_UTIL_ROOT}/hibf/misc/helper" "${CMAKE_CURRENT_BINARY_DIR}/hibf/misc/helper")

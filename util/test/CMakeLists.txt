# --------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
# --------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.21)

project (raptor_utility_test)

set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

get_filename_component (RAPTOR_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
option (RAPTOR_USE_CCACHE_IN_PARENT_PROJECT "" ON)
add_subdirectory ("${RAPTOR_ROOT_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/raptor")
target_compile_options (raptor_raptor INTERFACE "-pedantic" "-Wall" "-Wextra" "-Werror")

get_filename_component (RAPTOR_UTIL_ROOT "${CMAKE_CURRENT_LIST_DIR}/../dummy" DIRECTORY)

message (STATUS "Adding applications...")
add_subdirectory ("${RAPTOR_UTIL_ROOT}/applications/src" "${CMAKE_CURRENT_BINARY_DIR}/applications")
message (STATUS "Adding iScience...")
add_subdirectory ("${RAPTOR_UTIL_ROOT}/iScience" "${CMAKE_CURRENT_BINARY_DIR}/iScience")
message (STATUS "Adding thresholding...")
add_subdirectory ("${RAPTOR_UTIL_ROOT}/iScience/thresholding" "${CMAKE_CURRENT_BINARY_DIR}/thresholding")

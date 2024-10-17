# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25...3.30)

project (raptor_utility_test)

set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

get_filename_component (Raptor_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
include ("${Raptor_SOURCE_DIR}/cmake/configuration.cmake")
add_subdirectory ("${Raptor_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/raptor")
target_compile_options (raptor_interface INTERFACE "-pedantic" "-Wall" "-Wextra")

option (RAPTOR_WITH_WERROR "Report compiler warnings as errors." ON)

if (RAPTOR_WITH_WERROR)
    target_compile_options (raptor_interface INTERFACE "-Werror")
    message (STATUS "Building tests with -Werror.")
endif ()

get_filename_component (RAPTOR_UTIL_ROOT "${CMAKE_CURRENT_LIST_DIR}/../dummy" DIRECTORY)

message (STATUS "Adding applications...")
add_subdirectory ("${RAPTOR_UTIL_ROOT}/applications/src" "${CMAKE_CURRENT_BINARY_DIR}/applications")
message (STATUS "Adding iScience...")
add_subdirectory ("${RAPTOR_UTIL_ROOT}/iScience" "${CMAKE_CURRENT_BINARY_DIR}/iScience")
message (STATUS "Adding thresholding...")
add_subdirectory ("${RAPTOR_UTIL_ROOT}/iScience/thresholding" "${CMAKE_CURRENT_BINARY_DIR}/thresholding")

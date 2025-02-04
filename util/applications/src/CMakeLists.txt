# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25...3.30)

project (raptor_utility_common LANGUAGES CXX)

if (TARGET utility_common)
    return ()
endif ()

set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

get_filename_component (Raptor_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/../../.." ABSOLUTE)
include ("${Raptor_SOURCE_DIR}/cmake/configuration.cmake")
add_subdirectory ("${Raptor_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/raptor")

add_library (utility_common INTERFACE)
target_link_libraries ("utility_common" INTERFACE raptor::interface)
target_include_directories ("utility_common" INTERFACE "../include")

add_executable ("hibf_sizes" hibf_sizes.cpp)
target_link_libraries ("hibf_sizes" "utility_common")

add_executable ("simulate_sequence" simulate_sequence.cpp)
target_link_libraries ("simulate_sequence" "utility_common")

add_subdirectory (Genome_Biology)

add_subdirectory (iScience)

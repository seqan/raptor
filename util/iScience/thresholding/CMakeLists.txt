# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25)

# Define the application name.
project (minimizer_thresholds LANGUAGES CXX C)

# Make Release default build type
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE
         Release
         CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo" FORCE
    )
endif ()

# Specify the directories where to store the built archives, libraries and executables
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Dependencies
get_filename_component (Raptor_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/../../.." ABSOLUTE)
include ("${Raptor_SOURCE_DIR}/cmake/configuration.cmake")

add_library (common_thresholding INTERFACE)
target_link_libraries ("common_thresholding" INTERFACE seqan3::seqan3 sharg::sharg)
target_include_directories ("common_thresholding" INTERFACE include)

add_executable (query src/query.cpp)
target_link_libraries ("query" "common_thresholding")

add_executable (random_data src/random_data.cpp)
target_link_libraries ("random_data" "common_thresholding")

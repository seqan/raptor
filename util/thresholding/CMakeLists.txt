# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.15)

# Define the application name.
project (minimizer_thresholds CXX)

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

# Dependency: SeqAn3.
set (SEQAN3_SUBMODULES_DIR "../../lib")
find_package (SeqAn3 QUIET REQUIRED HINTS ../../lib/seqan3/build_system)

add_library (common_thresholding INTERFACE)
target_link_libraries ("common_thresholding" INTERFACE seqan3::seqan3)
target_include_directories ("common_thresholding" INTERFACE include)

add_executable (query src/query.cpp)
target_link_libraries ("query" "common_thresholding")

add_executable (random_data src/random_data.cpp)
target_link_libraries ("random_data" "common_thresholding")

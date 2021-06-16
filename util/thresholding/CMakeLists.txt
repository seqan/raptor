# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/eseiler/minimizer_thresholds/blob/master/LICENSE
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.8)

# Define the application name.
project (minimizer_thresholds CXX)

# Make Release default build type
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release CACHE STRING
         "Choose the type of build, options are: Debug Release RelWithDebInfo"
         FORCE)
endif ()

# Specify the directories where to store the built archives, libraries and executables
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Dependency: SeqAn3.
find_package (SeqAn3 QUIET REQUIRED HINTS ../../lib/seqan3/build_system)

add_library (common INTERFACE)
target_link_libraries ("common" INTERFACE seqan3::seqan3)
target_include_directories ("common" INTERFACE include)

add_executable (query src/query.cpp)
target_link_libraries ("query" "common")

add_executable (random_data src/random_data.cpp)
target_link_libraries ("random_data" "common")

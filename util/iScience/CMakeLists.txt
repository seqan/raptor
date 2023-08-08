# --------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
# --------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.18)

find_path (RAPTOR_CLONE_DIR
           NAMES build_system/raptor-config.cmake
           HINTS "${CMAKE_CURRENT_LIST_DIR}/../.."
)

list (APPEND CMAKE_MODULE_PATH "${RAPTOR_CLONE_DIR}/build_system")

include (raptor-config-version)

project (raptor_utility_common LANGUAGES CXX C)

find_package (Raptor REQUIRED HINTS ${RAPTOR_CLONE_DIR}/build_system)

include (${RAPTOR_CLONE_DIR}/test/cmake/raptor_require_ccache.cmake)
raptor_require_ccache ()

set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Define cmake configuration flags to configure and build external projects with the same flags as specified for
# this project.
set (SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "--no-warn-unused-cli")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_C_FLAGS=${CMAKE_C_FLAGS} -w")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS} -w")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}")

# Optionally build mason and DREAM-Yara.
option (RAPTOR_UTILITY_BUILD_MASON "Build mason. Needed to generate artificial data set." OFF)
option (RAPTOR_UTILITY_BUILD_DREAM_YARA "Build DREAM-Yara. Needed to map reads to the NCBI RefSeq." OFF)

if (RAPTOR_UTILITY_BUILD_MASON)
    message ("-- Configuring Mason")
    include (cmake/raptor_require_mason.cmake)
    raptor_require_mason ()
    message ("-- Configuring Mason - done")
else ()
    message ("-- Mason will not be built. Enable via 'cmake -DRAPTOR_UTILITY_BUILD_MASON=1 .'")
endif ()

if (RAPTOR_UTILITY_BUILD_DREAM_YARA)
    message ("-- Configuring DREAM-Yara")
    include (cmake/raptor_require_dream_yara.cmake)
    raptor_require_dream_yara ()
    message ("-- Configuring DREAM-Yara - done")
else ()
    message ("-- DREAM-Yara will not be built. Enable via 'cmake -DRAPTOR_UTILITY_BUILD_DREAM_YARA=1 .'")
endif ()

option (RAPTOR_UTILITY_BUILD_iScience_APPS "" ON)
if (RAPTOR_UTILITY_BUILD_iScience_APPS)
    option (RAPTOR_UTILITY_BUILD_Genome_Biology "" OFF)
    add_subdirectory ("${RAPTOR_CLONE_DIR}/util/applications/src" "${CMAKE_CURRENT_BINARY_DIR}/applications")
endif ()

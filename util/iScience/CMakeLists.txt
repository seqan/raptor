# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25...3.30)

set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

get_filename_component (Raptor_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
include ("${Raptor_SOURCE_DIR}/cmake/configuration.cmake")
add_subdirectory ("${Raptor_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/raptor")

# Define cmake configuration flags to configure and build external projects with the same flags as specified for
# this project.
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "--no-warn-unused-cli")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-Wno-dev")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_C_COMPILER_LAUNCHER=${CMAKE_C_COMPILER_LAUNCHER}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_C_FLAGS=${CMAKE_C_FLAGS} -w")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_COMPILER_LAUNCHER=${CMAKE_CXX_COMPILER_LAUNCHER}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS} -w")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_EXPORT_COMPILE_COMMANDS=${CMAKE_EXPORT_COMPILE_COMMANDS}")

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

add_subdirectory ("${Raptor_SOURCE_DIR}/util/applications/src" "${CMAKE_CURRENT_BINARY_DIR}/applications")

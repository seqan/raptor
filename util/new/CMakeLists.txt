# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.11)

# Define the application name.
project (raptor_utility)

# Make Release default build type
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE
         Release
         CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo" FORCE
    )
endif ()

get_filename_component (RAPTOR_ROOT_DIR "../.." ABSOLUTE)

# Specify the directories where to store the built archives, libraries and executables
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Messages
string (ASCII 27 Esc)
set (FontBold "${Esc}[1m")
set (FontReset "${Esc}[m")

# Dependency: SeqAn3.
set (SEQAN3_SUBMODULES_DIR "${RAPTOR_ROOT_DIR}/lib")
find_package (SeqAn3 QUIET REQUIRED HINTS ${RAPTOR_ROOT_DIR}/lib/seqan3/build_system)

# Allow to include CMake scripts from raptor.
list (APPEND CMAKE_MODULE_PATH "${RAPTOR_ROOT_DIR}/cmake/")

# Use ccache.
include (raptor_require_ccache)
raptor_require_ccache ()

add_library (common INTERFACE)
target_link_libraries ("common" INTERFACE seqan3::seqan3)
target_include_directories ("common" INTERFACE ${RAPTOR_ROOT_DIR}/include)
target_include_directories ("common" INTERFACE ${RAPTOR_ROOT_DIR}/lib/robin-hood-hashing/src/include)
target_compile_options ("common" INTERFACE "-pedantic" "-Wall" "-Wextra")

add_executable ("normalise_mantis_output" applications/normalise_mantis_output.cpp)
target_link_libraries ("normalise_mantis_output" "common")
install (TARGETS "normalise_mantis_output" DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

add_executable ("compare_mantis_raptor_output" applications/compare_mantis_raptor_output.cpp)
target_link_libraries ("compare_mantis_raptor_output" "common")
install (TARGETS "compare_mantis_raptor_output" DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

message (STATUS "${FontBold}You can run `make install` to build the application.${FontReset}")

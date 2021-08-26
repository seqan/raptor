# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https:#github.com/seqan/raptor/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.8)

# Define the application name and version.
project (raptor VERSION 2.0.0)

# Messages
string (ASCII 27 Esc)
set (FontBold "${Esc}[1m")
set (FontReset "${Esc}[m")

# Fallback to these values if there is no git or no git repository
set (RAPTOR_COMMIT_DATE "2021-08-20--no-git"
     CACHE STRING
     "Set to provide a commit date if git is not available or the source directory is not a git repository.")
set (RAPTOR_COMMIT_HASH "74f815358db47037e93a56b826a9df3692e55680--no-git"
     CACHE STRING
     "Set to provide a commit hash if git is not available or the source directory is not a git repository.")

# Extract git commit hash and date
find_package (Git QUIET)

if (GIT_FOUND)
    execute_process (COMMAND "${GIT_EXECUTABLE}" -C "${CMAKE_SOURCE_DIR}" rev-parse
                     WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
                     RESULT_VARIABLE is_no_git_repository
                     ERROR_QUIET)

    if (NOT is_no_git_repository)
        execute_process (COMMAND "${GIT_EXECUTABLE}" describe --always --abbrev=40 --dirty
                         WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
                         OUTPUT_VARIABLE RAPTOR_COMMIT_HASH
                         ERROR_QUIET
                         OUTPUT_STRIP_TRAILING_WHITESPACE)

        execute_process (COMMAND
                         "${GIT_EXECUTABLE}" log -1 --format=%ad --date=short
                         WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
                         OUTPUT_VARIABLE RAPTOR_COMMIT_DATE
                         ERROR_QUIET
                         OUTPUT_STRIP_TRAILING_WHITESPACE)
    endif ()
endif ()

## BUILD

# Make Release default build type
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release CACHE STRING
         "Choose the type of build, options are: Debug Release Coverage RelWithDebInfo MinSizeRel."
         FORCE)
endif ()

set (RAPTOR_NATIVE_BUILD ON CACHE BOOL "Optimize build for current architecture.")
if (RAPTOR_NATIVE_BUILD)
    message (STATUS "${FontBold}Native build enabled. Built binaries will be optimized for this system.${FontReset}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
else ()
     message (STATUS "${FontBold}Native build disabled. Detecting popcnt support.${FontReset}")
     include (CheckCXXCompilerFlag)
     check_cxx_compiler_flag ("-mpopcnt" RAPTOR_HAS_POPCNT)
     if (RAPTOR_HAS_POPCNT)
          set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mpopcnt")
     endif ()
endif ()

# Specify the directories where to store the built archives, libraries and executables
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Dependency: SeqAn3.
set (SEQAN3_CEREAL ON CACHE BOOL "Require cereal to be present.")
set (SEQAN3_SUBMODULES_DIR "lib")
find_package (SeqAn3 QUIET REQUIRED HINTS lib/seqan3/build_system)

# Use ccache.
set (RAPTOR_USE_CCACHE ON CACHE BOOL "Use ccache if available.")
if (RAPTOR_USE_CCACHE)
     include ("${SEQAN3_CLONE_DIR}/test/cmake/seqan3_require_ccache.cmake")
     seqan3_require_ccache ()
endif ()

# Add the application.
add_subdirectory (src)
message (STATUS "${FontBold}You can run `make` to build the application.${FontReset}")

## DOCUMENTATION

add_subdirectory (doc EXCLUDE_FROM_ALL)

## TEST

enable_testing ()
add_subdirectory (test EXCLUDE_FROM_ALL)

# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.11)

# Define the application name and version.
project (raptor VERSION 2.0.1)

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

set (RAPTOR_SUBMODULES_DIR "${CMAKE_CURRENT_LIST_DIR}/lib" CACHE STRING "Directory containing submodules.")

include (CheckCXXCompilerFlag)

set (RAPTOR_NATIVE_BUILD ON CACHE BOOL "Optimize build for current architecture.")
if ("${CMAKE_BUILD_TYPE}" MATCHES "Debug" OR "${CMAKE_BUILD_TYPE}" MATCHES "Coverage")
     message (STATUS "${FontBold}Native build disabled due to Debug/Coverage build.${FontReset}")
elseif (RAPTOR_NATIVE_BUILD)
     message (STATUS "${FontBold}Native build enabled. Built binaries will be optimized for this system.${FontReset}")
          set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
else ()
     message (STATUS "${FontBold}Native build disabled. Detecting popcnt support.${FontReset}")
     check_cxx_compiler_flag ("-mpopcnt" RAPTOR_HAS_POPCNT)
     if (RAPTOR_HAS_POPCNT)
          set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mpopcnt")
     endif ()
endif ()

check_cxx_compiler_flag ("-fopenmp" RAPTOR_HAS_OPENMP)
if (RAPTOR_HAS_OPENMP)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif ()

check_cxx_compiler_flag ("-fopenmp-simd" RAPTOR_HAS_OPENMP_SIMD)
if (RAPTOR_HAS_OPENMP_SIMD)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp-simd -DSIMDE_ENABLE_OPENMP")
endif ()

check_cxx_compiler_flag ("-fopenmp-simd" RAPTOR_SUPPRESS_GCC4_ABI)
if (RAPTOR_SUPPRESS_GCC4_ABI)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-psabi")
endif ()

# Specify the directories where to store the built archives, libraries and executables
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Dependency: SeqAn3.
set (SEQAN3_CEREAL ON CACHE BOOL "Require cereal to be present.")
set (SEQAN3_SUBMODULES_DIR "${RAPTOR_SUBMODULES_DIR}")
find_package (SeqAn3 QUIET REQUIRED HINTS ${RAPTOR_SUBMODULES_DIR}/seqan3/build_system)

# Define cmake configuration flags to configure and build external projects with the same flags as specified for
# this project.
set (SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_C_COMPILER_LAUNCHER=${CMAKE_C_COMPILER_LAUNCHER}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_COMPILER_LAUNCHER=${CMAKE_CXX_COMPILER_LAUNCHER}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}")

# Dependency: chopper
set (RAPTOR_CHOPPER_PATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libchopper_count_lib.a")
list (APPEND RAPTOR_CHOPPER_PATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libchopper_layout_lib.a")

if (EXISTS "${RAPTOR_CHOPPER_PATH}")
    add_custom_target (chopper_lib)
else ()
     set (RAPTOR_CHOPPER_ARGS "${SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS}")
     list (APPEND RAPTOR_CHOPPER_ARGS "-DCHOPPER_SUBMODULES_DIR=${RAPTOR_SUBMODULES_DIR}")
     list (APPEND RAPTOR_CHOPPER_ARGS "-DSEQAN3_SUBMODULES_DIR=${RAPTOR_SUBMODULES_DIR}")
     list (APPEND RAPTOR_CHOPPER_ARGS "-DCMAKE_INSTALL_LIBDIR=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
     list (APPEND RAPTOR_CHOPPER_ARGS "-DCMAKE_INSTALL_BINDIR=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

    include (ExternalProject)
    ExternalProject_Add (
        chopper_lib
        PREFIX chopper_lib
        SOURCE_DIR "${RAPTOR_SUBMODULES_DIR}/chopper"
        DOWNLOAD_COMMAND ""
        CMAKE_ARGS "${RAPTOR_CHOPPER_ARGS}"
        BUILD_BYPRODUCTS "${RAPTOR_CHOPPER_PATH}"
        UPDATE_DISCONNECTED 1
    )
endif ()

add_library (chopper_count STATIC IMPORTED)
add_dependencies (chopper_count chopper_lib)
set_target_properties (chopper_count PROPERTIES IMPORTED_LOCATION "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libchopper_count_lib.a")

add_library (chopper_layout STATIC IMPORTED)
add_dependencies (chopper_layout chopper_lib)
set_target_properties (chopper_layout PROPERTIES IMPORTED_LOCATION "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libchopper_layout_lib.a")

add_library (chopper INTERFACE IMPORTED)
set_property (TARGET chopper PROPERTY INTERFACE_LINK_LIBRARIES chopper_count chopper_layout)

# Allow to include CMake scripts from seqan3.
list (APPEND CMAKE_MODULE_PATH "${SEQAN3_CLONE_DIR}/test/cmake/")

# Use ccache.
set (RAPTOR_USE_CCACHE ON CACHE BOOL "Use ccache if available.")
if (RAPTOR_USE_CCACHE)
     include (seqan3_require_ccache)
     seqan3_require_ccache ()
endif ()

# Add the application.
add_subdirectory (src)
message (STATUS "${FontBold}You can run `make` to build the application.${FontReset}")

## DOCUMENTATION

add_subdirectory (doc EXCLUDE_FROM_ALL)

## TEST

set (RAPTOR_ENABLE_BENCHMARK OFF CACHE BOOL "Compile benchmarks as cli tests.")

enable_testing ()
add_subdirectory (test EXCLUDE_FROM_ALL)

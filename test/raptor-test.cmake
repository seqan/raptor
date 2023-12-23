# --------------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
# --------------------------------------------------------------------------------------------------------

# This file provides functionality common to the different test modules used by
# SeqAn3. To build tests, run cmake on one of the sub-folders in this directory
# which contain a CMakeLists.txt.

cmake_minimum_required (VERSION 3.10)

# require Raptor package
find_package (Raptor REQUIRED HINTS ${CMAKE_CURRENT_LIST_DIR}/../cmake)
include (${CMAKE_CURRENT_LIST_DIR}/../cmake/raptor-config-version.cmake)

include (CheckCXXSourceCompiles)
include (FindPackageHandleStandardArgs)
include (FindPackageMessage)

option (RAPTOR_TEST_BUILD_OFFLINE "Skip the update step of external projects." OFF)

message (STATUS "${ColourBold}Configuring tests${ColourReset}")

set (CPM_INDENT "  CMake Package Manager CPM: ")
CPMUsePackageLock ("${CMAKE_CURRENT_LIST_DIR}/../cmake/package-lock.cmake")

# ----------------------------------------------------------------------------
# Paths to folders.
# ----------------------------------------------------------------------------

find_path (RAPTOR_TEST_CMAKE_MODULE_DIR
           NAMES raptor_require_ccache.cmake
           HINTS "${CMAKE_CURRENT_LIST_DIR}/cmake/"
)
list (APPEND CMAKE_MODULE_PATH "${RAPTOR_TEST_CMAKE_MODULE_DIR}")

# ----------------------------------------------------------------------------
# Interface targets for the different test modules in seqan3.
# ----------------------------------------------------------------------------

enable_testing ()

# raptor::test exposes a base set of required flags, includes, definitions and
# libraries which are in common for **all** seqan3 tests
if (NOT TARGET raptor::test)
    add_library (raptor_test INTERFACE)
    target_compile_options (raptor_test INTERFACE "-pedantic" "-Wall" "-Wextra" "-Werror")

    # GCC12 and above: Disable warning about std::hardware_destructive_interference_size not being ABI-stable.
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 12)
            target_compile_options (raptor_test INTERFACE "-Wno-interference-size")
        endif ()
    endif ()

    target_link_libraries (raptor_test INTERFACE "raptor_lib")

    target_include_directories (raptor_test INTERFACE "${CMAKE_CURRENT_LIST_DIR}/include")

    add_library (raptor::test ALIAS raptor_test)
endif ()

# raptor::test::performance specifies required flags, includes and libraries
# needed for performance test cases in raptor/test/performance
if (NOT TARGET raptor::test::performance)
    add_library (raptor_test_performance INTERFACE)
    target_link_libraries (raptor_test_performance INTERFACE "raptor::test" "benchmark::benchmark_main")
    add_library (raptor::test::performance ALIAS raptor_test_performance)
endif ()

# raptor::test::unit specifies required flags, includes and libraries
# needed for unit test cases in raptor/test/unit
if (NOT TARGET raptor::test::unit)
    add_library (raptor_test_unit INTERFACE)

    # GCC12 has some bogus warnings. They will not be fixed in googletest.
    # https://github.com/google/googletest/issues/4232
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 12 AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 13)
            target_compile_options (raptor_test INTERFACE "-Wno-restrict")
        endif ()
    endif ()

    target_link_libraries (raptor_test_unit INTERFACE "raptor::test" "GTest::gtest_main")
    add_library (raptor::test::unit ALIAS raptor_test_unit)
endif ()

# raptor::test::header specifies required flags, includes and libraries
# needed for header test cases in raptor/test/header
if (NOT TARGET raptor::test::header)
    add_library (raptor_test_header INTERFACE)
    target_link_libraries (raptor_test_header INTERFACE "raptor::test::unit")
    target_link_libraries (raptor_test_header INTERFACE "raptor::test::performance")
    add_library (raptor::test::header ALIAS raptor_test_header)
endif ()

# ----------------------------------------------------------------------------
# Commonly used macros for the different test modules in seqan3.
# ----------------------------------------------------------------------------

include (app_datasources)
include (app_internal_datasources)
include (raptor_add_benchmark)
include (raptor_add_unit_test)
include (${CMAKE_CURRENT_LIST_DIR}/data/datasources.cmake)

# ----------------------------------------------------------------------------
# Add app.
# ----------------------------------------------------------------------------

get_filename_component (RAPTOR_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/../src" ABSOLUTE)
add_subdirectory ("${RAPTOR_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/src")

# ----------------------------------------------------------------------------
# Set directories for test output files, input data and binaries.
# ----------------------------------------------------------------------------

file (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output)
add_definitions (-DOUTPUTDIR=\"${CMAKE_CURRENT_BINARY_DIR}/output/\")
add_definitions (-DDATADIR=\"${CMAKE_CURRENT_BINARY_DIR}/data/\")
add_definitions (-DBINDIR=\"${CMAKE_BINARY_DIR}/bin/\")

# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25)

# ----------------------------------------------------------------------------
# Short-circuit if tests are already configured
# ----------------------------------------------------------------------------

if (TARGET raptor::test)
    return ()
endif ()

message (STATUS "${ColourBold}Configuring tests${ColourReset}")

# ----------------------------------------------------------------------------
# Add Raptor
# ----------------------------------------------------------------------------

get_filename_component (RAPTOR_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
add_subdirectory ("${RAPTOR_ROOT_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/raptor")
set_property (TARGET raptor PROPERTY RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
target_compile_options (raptor_raptor INTERFACE "-pedantic" "-Wall" "-Wextra" "-Werror")

# ----------------------------------------------------------------------------
# CPM
# ----------------------------------------------------------------------------

set (CPM_INDENT "CMake Package Manager CPM: ")
CPMUsePackageLock ("${CMAKE_CURRENT_LIST_DIR}/../cmake/package-lock.cmake")

# ----------------------------------------------------------------------------
# Paths to cmake modules
# ----------------------------------------------------------------------------

find_path (RAPTOR_TEST_CMAKE_MODULE_DIR
           NAMES app_datasources.cmake
           HINTS "${CMAKE_CURRENT_LIST_DIR}/cmake/"
)
list (APPEND CMAKE_MODULE_PATH "${RAPTOR_TEST_CMAKE_MODULE_DIR}")

# ----------------------------------------------------------------------------
# Interface targets for the different test modules
# ----------------------------------------------------------------------------

enable_testing ()

# ----------------------------------------------------------------------------
# raptor::test
# ----------------------------------------------------------------------------

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

# ----------------------------------------------------------------------------
# raptor::test::performance
# ----------------------------------------------------------------------------

add_library (raptor_test_performance INTERFACE)
target_link_libraries (raptor_test_performance INTERFACE "raptor::test" "benchmark::benchmark_main")
add_library (raptor::test::performance ALIAS raptor_test_performance)

# ----------------------------------------------------------------------------
# raptor::test::unit
# ----------------------------------------------------------------------------

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

# ----------------------------------------------------------------------------
# raptor::test::header
# ----------------------------------------------------------------------------

add_library (raptor_test_header INTERFACE)
target_link_libraries (raptor_test_header INTERFACE "raptor::test::unit")
target_link_libraries (raptor_test_header INTERFACE "raptor::test::performance")
add_library (raptor::test::header ALIAS raptor_test_header)

# ----------------------------------------------------------------------------
# Commonly used macros for the different test modules in seqan3.
# ----------------------------------------------------------------------------

include (app_datasources)
include (app_internal_datasources)
include (raptor_add_benchmark)
include (raptor_add_unit_test)
include (${CMAKE_CURRENT_LIST_DIR}/data/datasources.cmake)

# ----------------------------------------------------------------------------
# Set directories for test output files, input data and binaries.
# ----------------------------------------------------------------------------

file (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output)
add_definitions (-DOUTPUTDIR=\"${CMAKE_CURRENT_BINARY_DIR}/output/\")
add_definitions (-DDATADIR=\"${CMAKE_CURRENT_BINARY_DIR}/data/\")
add_definitions (-DBINDIR=\"${CMAKE_BINARY_DIR}/bin/\")

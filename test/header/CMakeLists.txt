# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25...3.30)

project (raptor_test_header LANGUAGES CXX C)

include (../raptor-test.cmake)
include (GoogleTest OPTIONAL)

CPMGetPackage (benchmark)
if (benchmark_ADDED)
    set_target_properties (benchmark PROPERTIES CXX_STANDARD ${CMAKE_CXX_STANDARD})
endif ()
CPMGetPackage (GTest)

find_path (SEQAN3_TEST_CMAKE_MODULE_DIR
           NAMES seqan3_test_files.cmake
           HINTS "${seqan3_SOURCE_DIR}/test/cmake/"
)
list (APPEND CMAKE_MODULE_PATH "${SEQAN3_TEST_CMAKE_MODULE_DIR}")

include (seqan3_test_files)
include (seqan3_test_component)

add_library ("raptor_header_test_lib" INTERFACE)
target_link_libraries ("raptor_header_test_lib" INTERFACE raptor::test raptor::test::header)
set (target "header_test")

set (header_base_path "${CMAKE_CURRENT_LIST_DIR}/../../include/")

file (GLOB_RECURSE header_files
      RELATIVE "${header_base_path}"
      "${header_base_path}/*.hpp"
)

list (FILTER header_files EXCLUDE REGEX
      "raptor/version.hpp|hibf/contrib/std|raptor/search/fpga/min_ibf_fpga_oneapi.hpp"
)

file (WRITE "${PROJECT_BINARY_DIR}/${target}.cpp" "")
add_executable (${target} ${PROJECT_BINARY_DIR}/${target}.cpp)
target_link_libraries (${target} raptor_header_test_lib)
add_test (NAME "header/${target}" COMMAND ${target})

foreach (header ${header_files})
    seqan3_test_component (header_test_name "${header}" TEST_NAME)
    seqan3_test_component (header_target_name "${header}" TARGET_UNIQUE_NAME)

    foreach (header_sub_test "header-guard" "no-self-include")
        set (header_target_source
             "${PROJECT_BINARY_DIR}/${target}_files/${header_test_name}.hpp-${header_sub_test}.cpp"
        )
        set (header_target "${target}--${header_target_name}-${header_sub_test}")

        string (REPLACE "-" "__" header_test_name_safe "${target}, ${header_target}")

        # we use add_custom_command to detect changes to a header file, which will update the generated source file
        add_custom_command (OUTPUT "${header_target_source}"
                            COMMAND "${CMAKE_COMMAND}" "-DHEADER_FILE_ABSOLUTE=${header_base_path}/${header}"
                                    "-DHEADER_FILE_INCLUDE=${header}" "-DHEADER_TARGET_SOURCE=${header_target_source}"
                                    "-DHEADER_TEST_NAME_SAFE=${header_test_name_safe}" "-DHEADER_COMPONENT=Raptor"
                                    "-DHEADER_SUB_TEST=${header_sub_test}" "-P"
                                    "${SEQAN3_CLONE_DIR}/test/header/generate_header_source.cmake"
                            DEPENDS "${header_base_path}/${header}"
                                    "${SEQAN3_CLONE_DIR}/test/header/generate_header_source.cmake"
        )

        add_library (${header_target} OBJECT "${header_target_source}")

        target_link_libraries (${header_target} raptor_header_test_lib)

        target_sources (${target} PRIVATE $<TARGET_OBJECTS:${header_target}>)
    endforeach ()
endforeach ()

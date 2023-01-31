# --------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
# --------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.10)

project (raptor_test_header CXX)

include (../raptor-test.cmake)
include (GoogleTest OPTIONAL)

raptor_require_benchmark ()
raptor_require_ccache ()
raptor_require_test ()

find_path (SEQAN3_TEST_CMAKE_MODULE_DIR
           NAMES seqan3_test_files.cmake
           HINTS "${RAPTOR_SUBMODULES_DIR}/seqan3/test/cmake/"
)
list (APPEND CMAKE_MODULE_PATH "${SEQAN3_TEST_CMAKE_MODULE_DIR}")

include (seqan3_test_files)
include (seqan3_test_component)

add_library ("raptor_header_test_lib" INTERFACE)
target_link_libraries ("raptor_header_test_lib" INTERFACE raptor::test::header)
set (target "header_test")

set (header_base_path "${CMAKE_CURRENT_LIST_DIR}/../../include/")

file (GLOB_RECURSE header_files
      RELATIVE "${header_base_path}"
      "${header_base_path}/*.hpp"
)

list (FILTER header_files EXCLUDE REGEX "raptor/version.hpp")

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

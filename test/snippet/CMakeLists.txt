# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25...3.30)
project (raptor_test_snippet LANGUAGES CXX C)

include (../raptor-test.cmake)
include (raptor_test_component)
include (raptor_test_files)

add_library (snippet_main snippet_main.cpp)
target_link_libraries (snippet_main PUBLIC "gtest_main" "gtest" "raptor::test")

macro (raptor_snippet test_name_prefix snippet snippet_base_path)
    raptor_test_component (snippet_target_name "${snippet}" TARGET_NAME)
    raptor_test_component (snippet_test_name "${snippet}" TEST_NAME)
    raptor_test_component (snippet_target_path "${snippet}" TARGET_PATH)
    set (target "${snippet_target_name}_snippet")

    add_executable (${target} "${snippet_base_path}/${snippet}")
    target_link_libraries (${target} PUBLIC snippet_main)
    set_target_properties (${target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                                "${CMAKE_CURRENT_BINARY_DIR}/${snippet_target_path}"
    )

    set (snippet_test_target "${test_name_prefix}/${snippet_test_name}_snippet")
    add_test (NAME "${snippet_test_target}" COMMAND ${target})

    # disable version checker, as it interferes with comparing the snippet output
    set_tests_properties ("${snippet_test_target}" PROPERTIES ENVIRONMENT SHARG_NO_VERSION_CHECK=0)

    set (snippet_compare_test_target "${test_name_prefix}/${snippet_test_name}_snippet_cmp_output")

    add_test (NAME "${snippet_compare_test_target}"
              COMMAND ${CMAKE_COMMAND} -DTARGET_FILE=$<TARGET_FILE:${target}>
                      -DSOURCE_FILE=${snippet_base_path}/${snippet} #
                      -P "${CMAKE_SOURCE_DIR}/compare_snippet_output.cmake"
    )

    # disable version checker, as it interferes with comparing the snippet output
    set_tests_properties ("${snippet_compare_test_target}" PROPERTIES ENVIRONMENT SHARG_NO_VERSION_CHECK=0)

    unset (snippet_target_name)
    unset (snippet_test_name)
    unset (snippet_target_path)
    unset (target)
endmacro ()

macro (raptor_snippets test_name_prefix snippet_base_path)
    raptor_test_files (snippet_files "${snippet_base_path}" "*.cpp")
    foreach (snippet ${snippet_files})
        raptor_snippet ("${test_name_prefix}" "${snippet}" "${snippet_base_path}")
    endforeach ()
endmacro ()

CPMGetPackage (GTest)

raptor_snippets ("snippet" "${CMAKE_CURRENT_SOURCE_DIR}")
raptor_snippets ("doc/snippet" "${CMAKE_CURRENT_SOURCE_DIR}/../../doc")

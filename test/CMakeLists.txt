cmake_minimum_required (VERSION 3.8)

# Set directories for test output files, input data and binaries.
file (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output)
add_definitions (-DOUTPUTDIR=\"${CMAKE_CURRENT_BINARY_DIR}/output/\")
add_definitions (-DDATADIR=\"${CMAKE_CURRENT_BINARY_DIR}/data/\")
add_definitions (-DBINDIR=\"${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/\")

# Define cmake configuration flags to configure and build external projects with the same flags as specified for
# this project.
set (SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}")
set (SEQAN3_TEST_CLONE_DIR "${PROJECT_BINARY_DIR}/vendor/googletest")

include ("${SEQAN3_CLONE_DIR}/test/cmake/seqan3_require_test.cmake")

seqan3_require_test ()

# Build tests just before their execution, because they have not been built with "all" target.
# The trick is here to provide a cmake file as a directory property that executes the build command.
file (WRITE "${CMAKE_CURRENT_BINARY_DIR}/build_test_targets.cmake"
            "execute_process(COMMAND ${CMAKE_COMMAND} --build . --target api_test)\n"
            "execute_process(COMMAND ${CMAKE_COMMAND} --build . --target cli_test)")
set_directory_properties (PROPERTIES TEST_INCLUDE_FILE "${CMAKE_CURRENT_BINARY_DIR}/build_test_targets.cmake")

# Define the test targets. All depending targets are built just before the test execution.
add_custom_target (api_test)
add_custom_target (cli_test)

# Test executables and libraries should not mix with the application files.
unset (CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
unset (CMAKE_LIBRARY_OUTPUT_DIRECTORY)
unset (CMAKE_RUNTIME_OUTPUT_DIRECTORY)

# A macro that adds an api or cli test.
macro (add_app_test test_filename test_alternative)
    # Extract the test target name.
    file (RELATIVE_PATH source_file "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_LIST_DIR}/${test_filename}")
    get_filename_component (target "${source_file}" NAME_WE)

    # Create the test target.
    add_executable (${target} ${test_filename})
    target_link_libraries (${target} "${PROJECT_NAME}_lib" seqan3::seqan3 gtest gtest_main)

    # Make seqan3::test available for both cli and api tests.
    target_include_directories(${target} PUBLIC "${SEQAN3_CLONE_DIR}/test/include")
    target_include_directories(${target} PUBLIC "${SEQAN3_TEST_CLONE_DIR}/googletest/include/")

    # Add the test to its general target (cli or api).
    if (${test_alternative} STREQUAL "CLI_TEST")
        add_dependencies (${target} "${PROJECT_NAME}") # cli test needs the application executable
        add_dependencies (cli_test ${target})
    elseif (${test_alternative} STREQUAL "API_TEST")
        add_dependencies (api_test ${target})
    endif ()

    # Generate and set the test name.
    get_filename_component (target_relative_path "${source_file}" DIRECTORY)
    if (target_relative_path)
        set (test_name "${target_relative_path}/${target}")
    else ()
        set (test_name "${target}")
    endif ()
    add_test (NAME "${test_name}" COMMAND ${target})

    unset (source_file)
    unset (target)
    unset (test_name)
endmacro ()

# A macro that adds an api test.
macro (add_api_test test_filename)
    add_app_test (${test_filename} API_TEST)
endmacro ()

# A macro that adds a cli test.
macro (add_cli_test test_filename)
    add_app_test (${test_filename} CLI_TEST)
endmacro ()

# Fetch data and add the tests.
include (data/datasources.cmake)
add_subdirectory (api)
add_subdirectory (cli)
add_subdirectory (coverage)

message (STATUS "${FontBold}You can run `make test` to build and run tests.${FontReset}")

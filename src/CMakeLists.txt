# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25...3.30)

if (TARGET raptor)
    return ()
endif ()

# Fallback to these values if there is no git or no git repository
set (RAPTOR_COMMIT_DATE
     "Unavailable"
     CACHE STRING
           "Set to provide a commit date if git is not available or the source directory is not a git repository."
)
set (RAPTOR_COMMIT_HASH
     "commit unavailable"
     CACHE STRING
           "Set to provide a commit hash if git is not available or the source directory is not a git repository."
)

# Extract git commit hash and date
find_package (Git QUIET)

if (GIT_FOUND)
    execute_process (COMMAND "${GIT_EXECUTABLE}" -C "${CMAKE_CURRENT_LIST_DIR}" rev-parse
                     WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}"
                     RESULT_VARIABLE is_no_git_repository
                     ERROR_QUIET
    )

    if (NOT is_no_git_repository)
        execute_process (COMMAND "${GIT_EXECUTABLE}" describe --always --abbrev=40 --dirty
                         WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}"
                         OUTPUT_VARIABLE RAPTOR_COMMIT_HASH
                         ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE
        )

        execute_process (COMMAND "${GIT_EXECUTABLE}" log -1 --format=%ad --date=short
                         WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}"
                         OUTPUT_VARIABLE RAPTOR_COMMIT_DATE
                         ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    endif ()
endif ()

# Raptor library
add_library ("raptor_lib" INTERFACE)
target_link_libraries ("raptor_lib"
                       INTERFACE "raptor::argument_parsing"
                                 "raptor::build"
                                 "raptor::prepare"
                                 "raptor::search"
                                 "raptor::threshold"
                                 "raptor::upgrade"
                                 "raptor::layout"
                                 "raptor::update"
)
add_library (raptor::raptor ALIAS raptor_lib)

# Raptor executable
add_executable ("raptor" raptor.cpp)
target_link_libraries ("raptor" PUBLIC "raptor::raptor")
set_property (TARGET raptor PROPERTY RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    target_compile_options (raptor PUBLIC "-fwhole-program")
endif ()
set_source_files_properties (raptor.cpp
                             PROPERTIES COMPILE_DEFINITIONS
                                        "RAPTOR_VERSION=\"${PACKAGE_VERSION} (${RAPTOR_COMMIT_HASH})\";RAPTOR_DATE=\"${RAPTOR_COMMIT_DATE}\""
)

add_subdirectory (argument_parsing)
add_subdirectory (build)
add_subdirectory (layout)
add_subdirectory (search)
add_subdirectory (prepare)
add_subdirectory (threshold)
add_subdirectory (upgrade)
add_subdirectory (update)

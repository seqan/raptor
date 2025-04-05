# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25...3.30)

if (TARGET raptor::argument_parsing)
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

add_library ("raptor_argument_parsing" STATIC
             build_arguments.cpp
             build_parsing.cpp
             compute_bin_size.cpp
             init_shared_meta.cpp
             parse_bin_path.cpp
             prepare_parsing.cpp
             search_arguments.cpp
             search_parsing.cpp
             update_parsing.cpp
             upgrade_parsing.cpp
)

set_source_files_properties (search_parsing.cpp
                             PROPERTIES COMPILE_DEFINITIONS
                                        "FPGA_WINDOWS=${WINDOW_SIZE_STRING};FPGA_KMERS=${MIN_IBF_K_STRING};FPGA_BINS=${BIN_COUNT_STRING};FPGA_KERNELS=${KERNEL_COPYS_STRING}"
)

set_source_files_properties (init_shared_meta.cpp
                             PROPERTIES COMPILE_DEFINITIONS
                                        "RAPTOR_VERSION=\"${PACKAGE_VERSION} (${RAPTOR_COMMIT_HASH})\";RAPTOR_DATE=\"${RAPTOR_COMMIT_DATE}\""
)

target_link_libraries ("raptor_argument_parsing" PUBLIC "raptor::interface")
add_library (raptor::argument_parsing ALIAS raptor_argument_parsing)

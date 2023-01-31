cmake_minimum_required (VERSION 3.18)

if (NOT TARGET raptor_argument_parsing)
    # Fallback to these values if there is no git or no git repository
    set (RAPTOR_COMMIT_DATE
         "2022-01-17--no-git"
         CACHE STRING
               "Set to provide a commit date if git is not available or the source directory is not a git repository."
    )
    set (RAPTOR_COMMIT_HASH
         "80736310204cb2fa5724f1d9316e6391667dfe08--no-git"
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
                 ../build/hibf/parse_chopper_pack_line.cpp
                 build_parsing.cpp
                 compute_bin_size.cpp
                 init_shared_meta.cpp
                 parse_bin_path.cpp
                 prepare_parsing.cpp
                 search_parsing.cpp
                 upgrade_parsing.cpp
    )

    target_compile_definitions ("raptor_argument_parsing"
                                PUBLIC "-DRAPTOR_VERSION=\"${RAPTOR_PROJECT_VERSION} (${RAPTOR_COMMIT_HASH})\""
                                       "-DRAPTOR_DATE=\"${RAPTOR_COMMIT_DATE}\""
    )

    target_link_libraries ("raptor_argument_parsing" PUBLIC "raptor_interface")
endif ()

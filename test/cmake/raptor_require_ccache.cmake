# --------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
# --------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.18)

include (FindPackageMessage)

# Uses `ccache` to cache build results.
#
# See also
# * https://ccache.dev/
# * https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_LAUNCHER.html
macro (raptor_require_ccache)
    set (RAPTOR_USE_CCACHE
         ON
         CACHE BOOL "Use ccache if available."
    )
    set (RAPTOR_FPROFILE_ABS_PATH "-fprofile-abs-path")
    if (RAPTOR_USE_CCACHE)
        find_program (CCACHE_PROGRAM ccache)

        if (NOT CCACHE_PROGRAM)
            find_package_message (CCACHE_PROGRAM "  Finding program ccache - Failed" "[${CCACHE_PROGRAM}]")
        else ()
            find_package_message (CCACHE_PROGRAM "  Finding program ccache - Success" "[${CCACHE_PROGRAM}]")
            set (RAPTOR_FPROFILE_ABS_PATH "--ccache-skip -fprofile-abs-path")
            # New option since cmake >= 3.4:
            # https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_LAUNCHER.html
            if (NOT CMAKE_VERSION VERSION_LESS 3.15) # cmake >= 3.15
                list (PREPEND CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
                list (PREPEND CMAKE_C_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
            else ()
                # prepend ccache to CMAKE_CXX_COMPILER_LAUNCHER
                list (INSERT CMAKE_CXX_COMPILER_LAUNCHER 0 "${CCACHE_PROGRAM}")
                list (INSERT CMAKE_C_COMPILER_LAUNCHER 0 "${CCACHE_PROGRAM}")
            endif ()

            # use ccache in external cmake projects
            list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS
                  "-DCMAKE_CXX_COMPILER_LAUNCHER=${CMAKE_CXX_COMPILER_LAUNCHER}"
            )
            list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_C_COMPILER_LAUNCHER=${CMAKE_C_COMPILER_LAUNCHER}")

            if (NOT CMAKE_VERSION VERSION_LESS 3.21) # cmake >= 3.21
                list (PREPEND CMAKE_CXX_LINKER_LAUNCHER "${CCACHE_PROGRAM}")
                list (PREPEND CMAKE_C_LINKER_LAUNCHER "${CCACHE_PROGRAM}")
                list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS
                      "-DCMAKE_CXX_LINKER_LAUNCHER=${CMAKE_CXX_LINKER_LAUNCHER}"
                )
                list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_C_LINKER_LAUNCHER=${CMAKE_C_LINKER_LAUNCHER}")
            else ()
                set_property (GLOBAL PROPERTY RULE_LAUNCH_LINK "${CCACHE_PROGRAM}")
            endif ()
        endif ()
        unset (CCACHE_PROGRAM)
    endif ()
endmacro ()

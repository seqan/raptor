# --------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
# --------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.21)

include (FindPackageMessage)

# Uses `ccache` to cache build results.
#
# See also
# * https://ccache.dev/
# * https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_LAUNCHER.html
macro (raptor_require_ccache)
    option (RAPTOR_USE_CCACHE "Use ccache if available." ON)
    option (RAPTOR_USE_CCACHE_IN_PARENT_PROJECT "Use ccache in parent project if available." OFF)
    set (RAPTOR_FPROFILE_ABS_PATH "-fprofile-abs-path")
    if (RAPTOR_USE_CCACHE)
        find_program (CCACHE_PROGRAM ccache)

        if (NOT CCACHE_PROGRAM)
            find_package_message (CCACHE_PROGRAM "  Ccache program:             not available" "[${CCACHE_PROGRAM}]")
        else ()
            find_package_message (CCACHE_PROGRAM "  Ccache program:             available" "[${CCACHE_PROGRAM}]")
            set (RAPTOR_FPROFILE_ABS_PATH "--ccache-skip -fprofile-abs-path")

            list (PREPEND CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
            list (PREPEND CMAKE_C_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")

            # use ccache in external cmake projects
            list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS
                  "-DCMAKE_CXX_COMPILER_LAUNCHER=${CMAKE_CXX_COMPILER_LAUNCHER}"
            )
            list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_C_COMPILER_LAUNCHER=${CMAKE_C_COMPILER_LAUNCHER}")

            list (PREPEND CMAKE_CXX_LINKER_LAUNCHER "${CCACHE_PROGRAM}")
            list (PREPEND CMAKE_C_LINKER_LAUNCHER "${CCACHE_PROGRAM}")
            list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_LINKER_LAUNCHER=${CMAKE_CXX_LINKER_LAUNCHER}")
            list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_C_LINKER_LAUNCHER=${CMAKE_C_LINKER_LAUNCHER}")
        endif ()

        if (RAPTOR_USE_CCACHE_IN_PARENT_PROJECT)
            set (RAPTOR_FPROFILE_ABS_PATH
                 ${RAPTOR_FPROFILE_ABS_PATH}
                 PARENT_SCOPE
            )
            set (SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS
                 ${SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS}
                 PARENT_SCOPE
            )
            set (CMAKE_CXX_COMPILER_LAUNCHER
                 ${CMAKE_CXX_COMPILER_LAUNCHER}
                 PARENT_SCOPE
            )
            set (CMAKE_C_COMPILER_LAUNCHER
                 ${CMAKE_C_COMPILER_LAUNCHER}
                 PARENT_SCOPE
            )
            set (CMAKE_CXX_LINKER_LAUNCHER
                 ${CMAKE_CXX_LINKER_LAUNCHER}
                 PARENT_SCOPE
            )
            set (CMAKE_C_LINKER_LAUNCHER
                 ${CMAKE_C_LINKER_LAUNCHER}
                 PARENT_SCOPE
            )
        endif ()
        unset (CCACHE_PROGRAM)
    endif ()
endmacro ()

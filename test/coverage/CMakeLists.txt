# --------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
# --------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.15)
cmake_minimum_required (VERSION 3.10)

project (raptor_test_coverage CXX)

include (../cmake/raptor_require_ccache.cmake)
raptor_require_ccache ()

# Add a custom build type: Coverage

set (CMAKE_CXX_FLAGS_COVERAGE
     "${CMAKE_CXX_FLAGS_DEBUG} --coverage -fprofile-arcs -ftest-coverage ${RAPTOR_FPROFILE_ABS_PATH}"
     CACHE STRING "Flags used by the C++ compiler during coverage builds." FORCE
)
set (CMAKE_C_FLAGS_COVERAGE
     "${CMAKE_C_FLAGS_DEBUG} --coverage -fprofile-arcs -ftest-coverage ${RAPTOR_FPROFILE_ABS_PATH}"
     CACHE STRING "Flags used by the C compiler during coverage builds." FORCE
)
set (CMAKE_EXE_LINKER_FLAGS_COVERAGE
     "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -Wl,-lgcov"
     CACHE STRING "Flags used for linking binaries during coverage builds." FORCE
)
set (CMAKE_SHARED_LINKER_FLAGS_COVERAGE
     "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} -Wl,-lgcov"
     CACHE STRING "Flags used by the shared libraries linker during coverage builds." FORCE
)

mark_as_advanced (CMAKE_CXX_FLAGS_COVERAGE CMAKE_C_FLAGS_COVERAGE CMAKE_EXE_LINKER_FLAGS_COVERAGE
                  CMAKE_SHARED_LINKER_FLAGS_COVERAGE
)

enable_testing ()
add_subdirectory ("../unit" "${CMAKE_CURRENT_BINARY_DIR}/unit")

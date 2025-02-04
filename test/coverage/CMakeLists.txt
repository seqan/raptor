# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25...3.30)

project (raptor_test_coverage LANGUAGES CXX C)

# Add a custom build type: Coverage

set (CMAKE_CXX_FLAGS_COVERAGE
     "${CMAKE_CXX_FLAGS_DEBUG} --coverage -fprofile-arcs -ftest-coverage -fprofile-abs-path"
     CACHE STRING "Flags used by the C++ compiler during coverage builds." FORCE
)
set (CMAKE_C_FLAGS_COVERAGE
     "${CMAKE_C_FLAGS_DEBUG} --coverage -fprofile-arcs -ftest-coverage -fprofile-abs-path"
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

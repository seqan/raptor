# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25...3.30)

if (TARGET raptor)
    return ()
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

add_subdirectory (argument_parsing)
add_subdirectory (build)
add_subdirectory (layout)
add_subdirectory (search)
add_subdirectory (prepare)
add_subdirectory (threshold)
add_subdirectory (upgrade)
add_subdirectory (update)

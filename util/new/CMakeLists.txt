# --------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
# --------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.18)

# Define the application name.
project (raptor_utility)

# Make Release default build type
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE
         Release
         CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo" FORCE
    )
endif ()

include (CheckCXXCompilerFlag)

set (RAPTOR_NATIVE_BUILD
     ON
     CACHE BOOL "Optimize build for current architecture."
)
if ("${CMAKE_BUILD_TYPE}" MATCHES "Debug" OR "${CMAKE_BUILD_TYPE}" MATCHES "Coverage")
    message (STATUS "${FontBold}Native build disabled due to Debug/Coverage build.${FontReset}")
elseif (RAPTOR_NATIVE_BUILD)
    message (STATUS "${FontBold}Native build enabled. Built binaries will be optimized for this system.${FontReset}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
else ()
    message (STATUS "${FontBold}Native build disabled. Detecting popcnt support.${FontReset}")
    check_cxx_compiler_flag ("-mpopcnt" RAPTOR_HAS_POPCNT)
    if (RAPTOR_HAS_POPCNT)
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mpopcnt")
    endif ()
endif ()

get_filename_component (RAPTOR_ROOT_DIR "../.." ABSOLUTE)

# Specify the directories where to store the built archives, libraries and executables
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Messages
string (ASCII 27 Esc)
set (FontBold "${Esc}[1m")
set (FontReset "${Esc}[m")

# Dependency: SeqAn3.
set (SEQAN3_SUBMODULES_DIR "${RAPTOR_ROOT_DIR}/lib")
find_package (SeqAn3 QUIET REQUIRED HINTS ${RAPTOR_ROOT_DIR}/lib/seqan3/build_system)
find_package (Sharg QUIET REQUIRED HINTS ${RAPTOR_ROOT_DIR}/lib/sharg-parser/build_system)

# Allow to include CMake scripts from raptor.
list (APPEND CMAKE_MODULE_PATH "${RAPTOR_ROOT_DIR}/test/cmake/")

# Use ccache.
include (raptor_require_ccache)
raptor_require_ccache ()

add_library (new_common INTERFACE)
target_link_libraries ("new_common" INTERFACE seqan3::seqan3 sharg::sharg)
target_include_directories ("new_common" INTERFACE ${RAPTOR_ROOT_DIR}/include)
target_include_directories ("new_common" INTERFACE ${RAPTOR_ROOT_DIR}/lib/robin-hood-hashing/src/include)
target_compile_options ("new_common" INTERFACE "-pedantic" "-Wall" "-Wextra")

add_executable ("normalise_mantis_output" applications/normalise_mantis_output.cpp)
target_link_libraries ("normalise_mantis_output" "new_common")

add_executable ("normalise_raptor_output" applications/normalise_raptor_output.cpp)
target_link_libraries ("normalise_raptor_output" "common")
install (TARGETS "normalise_raptor_output" DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

add_executable ("compare_mantis_raptor_output" applications/compare_mantis_raptor_output.cpp)
target_link_libraries ("compare_mantis_raptor_output" "new_common")

add_executable ("normalise_raptor_output" applications/normalise_raptor_output.cpp)
target_link_libraries ("normalise_raptor_output" "new_common")

message (STATUS "${FontBold}You can run `make install` to build the application.${FontReset}")

# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.11)

# Define the application name.
project (raptor_utility)

# Make Release default build type
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE
         Release
         CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo" FORCE
    )
endif ()

# Specify the directories where to store the built archives, libraries and executables
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Messages
string (ASCII 27 Esc)
set (FontBold "${Esc}[1m")
set (FontReset "${Esc}[m")

# Dependency: SeqAn3.
set (SEQAN3_SUBMODULES_DIR "../lib")
find_package (SeqAn3 QUIET REQUIRED HINTS ../lib/seqan3/build_system)

# Allow to include CMake scripts from seqan3.
list (APPEND CMAKE_MODULE_PATH "${SEQAN3_CLONE_DIR}/test/cmake/")

# Define cmake configuration flags to configure and build external projects with the same flags as specified for
# this project.
set (SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "--no-warn-unused-cli")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}")

# Use ccache.
include (seqan3_require_ccache)
seqan3_require_ccache ()

# Optionally build mason and DREAM-Yara.
option (RAPTOR_UTILITY_BUILD_MASON "Build mason. Needed to generate artificial data set." OFF)
option (RAPTOR_UTILITY_BUILD_DREAM_YARA "Build DREAM-Yara. Needed to map reads to the NCBI RefSeq." OFF)

if (RAPTOR_UTILITY_BUILD_MASON)
    message ("-- Configuring Mason")
    include (cmake/raptor_require_mason.cmake)
    raptor_require_mason ()
    message ("-- Configuring Mason - done")
else ()
    message ("-- Mason will not be built. Enable via 'cmake -DRAPTOR_UTILITY_BUILD_MASON=1 .'")
endif ()

if (RAPTOR_UTILITY_BUILD_DREAM_YARA)
    message ("-- Configuring DREAM-Yara")
    include (cmake/raptor_require_dream_yara.cmake)
    raptor_require_dream_yara ()
    message ("-- Configuring DREAM-Yara - done")
else ()
    message ("-- DREAM-Yara will not be built. Enable via 'cmake -DRAPTOR_UTILITY_BUILD_DREAM_YARA=1 .'")
endif ()

add_library (common INTERFACE)
target_link_libraries ("common" INTERFACE seqan3::seqan3)
target_include_directories ("common" INTERFACE ../include)
target_include_directories ("common" INTERFACE ../lib/robin-hood-hashing/src/include)
target_compile_options ("common" INTERFACE "-pedantic" "-Wall" "-Wextra")

add_executable ("generate_reads" src/applications/generate_reads.cpp)
target_link_libraries ("generate_reads" "common")
install (TARGETS "generate_reads" DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

add_executable ("split_sequence" src/applications/split_sequence.cpp)
target_link_libraries ("split_sequence" "common")
install (TARGETS "split_sequence" DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

add_executable ("count_minimiser" ../src/argument_parsing/parse_bin_path.cpp
                                  ../src/build/hibf/parse_chopper_pack_line.cpp src/applications/count_minimiser.cpp
)
target_link_libraries ("count_minimiser" "common")
install (TARGETS "count_minimiser" DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

add_executable ("threshold_info"
                ../src/threshold/multiple_error_model.cpp
                ../src/threshold/one_error_model.cpp
                ../src/threshold/one_indirect_error_model.cpp
                ../src/threshold/pascal_row.cpp
                ../src/threshold/precompute_correction.cpp
                ../src/threshold/precompute_threshold.cpp
                src/applications/threshold_info.cpp
)
target_link_libraries ("threshold_info" "common")
install (TARGETS "threshold_info" DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

add_executable ("apply_taxsbp" src/applications/apply_taxsbp.cpp)
target_link_libraries ("apply_taxsbp" "common")
install (TARGETS "apply_taxsbp" DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

add_executable ("generate_reads_refseq" src/applications/generate_reads_refseq.cpp)
target_link_libraries ("generate_reads_refseq" "common")
install (TARGETS "generate_reads_refseq" DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

message (STATUS "${FontBold}You can run `make install` to build the application.${FontReset}")

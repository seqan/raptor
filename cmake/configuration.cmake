# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25...3.30)

# ----------------------------------------------------------------------------
# Short-circuit if Raptor is already configured
# ----------------------------------------------------------------------------

if (TARGET raptor::interface)
    return ()
endif ()

# ----------------------------------------------------------------------------
# Greeter
# ----------------------------------------------------------------------------

message (STATUS "Finding Raptor (${RAPTOR_VERSION}) and checking requirements")

# ----------------------------------------------------------------------------
# Pretty printing and error handling
# ----------------------------------------------------------------------------

macro (raptor_config_print text)
    message (STATUS "  ${text}")
endmacro ()

macro (raptor_config_error text)
    message (FATAL_ERROR "  ${text}")
endmacro ()

# ----------------------------------------------------------------------------
# C++ standard
# ----------------------------------------------------------------------------

if (NOT DEFINED CMAKE_CXX_STANDARD)
    set (CMAKE_CXX_STANDARD 23)
endif ()

if (NOT DEFINED CMAKE_CXX_STANDARD_REQUIRED)
    set (CMAKE_CXX_STANDARD_REQUIRED OFF)
endif ()

if (NOT DEFINED CMAKE_CXX_EXTENSIONS)
    set (CMAKE_CXX_EXTENSIONS OFF)
endif ()

# ----------------------------------------------------------------------------
# LTO
# ----------------------------------------------------------------------------

include (CheckIPOSupported)
check_ipo_supported (RESULT RAPTOR_HAS_LTO OUTPUT RAPTOR_HAS_LTO_OUTPUT)
if (RAPTOR_HAS_LTO)
    set (CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif ()

# ----------------------------------------------------------------------------
# CPM
# ----------------------------------------------------------------------------

set (CPM_INDENT "  CMake Package Manager CPM: ")
include (${Raptor_SOURCE_DIR}/cmake/CPM.cmake)
CPMUsePackageLock (${Raptor_SOURCE_DIR}/cmake/package-lock.cmake)

# ----------------------------------------------------------------------------
# Find or add dependencies
# ----------------------------------------------------------------------------

CPMGetPackage (use_ccache)
CPMGetPackage (hibf)
CPMGetPackage (sharg)
CPMGetPackage (seqan3)
CPMGetPackage (chopper)

# ----------------------------------------------------------------------------
# Find Raptor include path
# ----------------------------------------------------------------------------

find_path (RAPTOR_INCLUDE_DIR
           NAMES raptor/version.hpp
           HINTS "${Raptor_SOURCE_DIR}/include"
)

if (RAPTOR_INCLUDE_DIR)
    raptor_config_print ("Raptor include dir found:   ${RAPTOR_INCLUDE_DIR}")
else ()
    raptor_config_error ("Raptor include directory could not be found (RAPTOR_INCLUDE_DIR: '${RAPTOR_INCLUDE_DIR}')")
endif ()

# ----------------------------------------------------------------------------
# Export targets
# ----------------------------------------------------------------------------

add_library (raptor_interface INTERFACE)
target_link_libraries (raptor_interface INTERFACE sharg::sharg seqan3::seqan3 seqan::hibf)
target_include_directories (raptor_interface INTERFACE "${RAPTOR_INCLUDE_DIR}")

get_target_property (CHOPPER_INCLUDE_DIR chopper::interface INTERFACE_INCLUDE_DIRECTORIES)
target_include_directories (raptor_interface SYSTEM INTERFACE "${CHOPPER_INCLUDE_DIR}")

# !Workaround: Get seqan3 test include dir from seqan3 target
find_path (SEQAN3_TEST_INCLUDE_DIR
           NAMES seqan3/test/tmp_directory.hpp
           HINTS "${seqan3_SOURCE_DIR}/test/include"
)
target_include_directories (raptor_interface SYSTEM INTERFACE "${SEQAN3_TEST_INCLUDE_DIR}")

add_library (raptor::interface ALIAS raptor_interface)

# ----------------------------------------------------------------------------
# FPGA
# ----------------------------------------------------------------------------

include (${Raptor_SOURCE_DIR}/cmake/fpga_config.cmake)

# SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25)

# ----------------------------------------------------------------------------
# Short-circuit if Raptor is already configured
# ----------------------------------------------------------------------------

if (TARGET raptor::raptor)
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
# CPM
# ----------------------------------------------------------------------------

set (CPM_INDENT "CMake Package Manager CPM: ")
include (CPM)
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

add_library (raptor_raptor INTERFACE)
target_link_libraries (raptor_raptor INTERFACE seqan::hibf sharg::sharg seqan3::seqan3)
target_include_directories (raptor_raptor INTERFACE "${RAPTOR_INCLUDE_DIR}")

# !Workaround: Get chopper include dir from chopper_shared target
find_path (CHOPPER_INCLUDE_DIR
           NAMES chopper/configuration.hpp
           HINTS "${chopper_SOURCE_DIR}/include"
)
target_include_directories (raptor_raptor SYSTEM INTERFACE "${CHOPPER_INCLUDE_DIR}")

# !Workaround: Get seqan3 test include dir from seqan3 target
find_path (SEQAN3_TEST_INCLUDE_DIR
           NAMES seqan3/test/tmp_directory.hpp
           HINTS "${seqan3_SOURCE_DIR}/test/include"
)
target_include_directories (raptor_raptor SYSTEM INTERFACE "${SEQAN3_TEST_INCLUDE_DIR}")

add_library (raptor::raptor ALIAS raptor_raptor)

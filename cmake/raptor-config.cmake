# --------------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
# --------------------------------------------------------------------------------------------------------
#
# This CMake module will try to find RAPTOR and its dependencies.  You can use
# it the same way you would use any other CMake module.
#
#   find_package (RAPTOR [REQUIRED] ...)
#
# Since this makes a difference for CMAKE, pay attention to the case
# ("RAPTOR", "RAPTOR" and "raptor" are all valid, but other names not).
#
# RAPTOR has the following platform requirements:
#
#   C++20
#   pthread
#
# RAPTOR has the following optional dependencies:
#
#   ZLIB      -- zlib compression library
#   BZip2     -- libbz2 compression library
#
# If you don't wish for these to be detected (and used), you may define RAPTOR_NO_ZLIB,
# RAPTOR_NO_BZIP2, and RAPTOR_NO_CEREAL, respectively.
#
# If you wish to require the presence of ZLIB or BZip2, just check for the module before
# finding RAPTOR, e.g. "find_package (ZLIB REQUIRED)".
# If you wish to require the presence of CEREAL, you may define RAPTOR_CEREAL.
#
# Once the search has been performed, the following variables will be set.
#
#   RAPTOR_FOUND            -- Indicate whether RAPTOR was found and requirements met.
#
#   RAPTOR_VERSION          -- The version as string, e.g. "3.0.0"
#   RAPTOR_VERSION_MAJOR    -- e.g. 3
#   RAPTOR_VERSION_MINOR    -- e.g. 0
#   RAPTOR_VERSION_PATCH    -- e.g. 0
#
#   RAPTOR_INCLUDE_DIRS     -- to be passed to include_directories ()
#   RAPTOR_LIBRARIES        -- to be passed to target_link_libraries ()
#   RAPTOR_DEFINITIONS      -- to be passed to add_definitions ()
#   RAPTOR_CXX_FLAGS        -- to be added to CMAKE_CXX_FLAGS
#
# Additionally, the following [IMPORTED][IMPORTED] targets are defined:
#
#   raptor::raptor          -- interface target where
#                                  target_link_libraries(target raptor::raptor)
#                              automatically sets
#                                  target_include_directories(target $RAPTOR_INCLUDE_DIRS),
#                                  target_link_libraries(target $RAPTOR_LIBRARIES),
#                                  target_compile_definitions(target $RAPTOR_DEFINITIONS) and
#                                  target_compile_options(target $RAPTOR_CXX_FLAGS)
#                              for a target.
#
#   [IMPORTED]: https://cmake.org/cmake/help/v3.10/prop_tgt/IMPORTED.html#prop_tgt:IMPORTED
#
# ============================================================================

cmake_minimum_required (VERSION 3.5...3.12)

# ----------------------------------------------------------------------------
# Greeter
# ----------------------------------------------------------------------------

message (STATUS "Finding Raptor and checking requirements")

# ----------------------------------------------------------------------------
# Includes
# ----------------------------------------------------------------------------

include (${RAPTOR_CLONE_DIR}/test/cmake/raptor_require_ccache.cmake)
raptor_require_ccache ()

set (CPM_INDENT "  CMake Package Manager CPM: ")
include (${CMAKE_CURRENT_LIST_DIR}/CPM.cmake)
CPMUsePackageLock (${CMAKE_CURRENT_LIST_DIR}/package-lock.cmake)

include (CheckIncludeFileCXX)
include (CheckCXXSourceCompiles)
include (FindPackageHandleStandardArgs)
include (CheckCXXCompilerFlag)

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
# Find or add dependencies
# ----------------------------------------------------------------------------

CPMGetPackage (hibf)
CPMGetPackage (sharg)
CPMGetPackage (seqan3)
CPMGetPackage (chopper)

# ----------------------------------------------------------------------------
# Find RAPTOR include path
# ----------------------------------------------------------------------------

# Note that raptor-config.cmake can be standalone and thus RAPTOR_CLONE_DIR might be empty.
# * `RAPTOR_CLONE_DIR` was already found in raptor-config-version.cmake
# * `RAPTOR_INCLUDE_DIR` was already found in raptor-config-version.cmake
find_path (RAPTOR_SUBMODULES_DIR
           NAMES seqan3
           HINTS "${RAPTOR_CLONE_DIR}/lib" "${RAPTOR_INCLUDE_DIR}/raptor"
)

if (RAPTOR_INCLUDE_DIR)
    raptor_config_print ("RAPTOR include dir found:   ${RAPTOR_INCLUDE_DIR}")
else ()
    raptor_config_error ("RAPTOR include directory could not be found (RAPTOR_INCLUDE_DIR: '${RAPTOR_INCLUDE_DIR}')")
endif ()

# ----------------------------------------------------------------------------
# Finish find_package call
# ----------------------------------------------------------------------------

find_package_handle_standard_args (${CMAKE_FIND_PACKAGE_NAME} REQUIRED_VARS RAPTOR_INCLUDE_DIR)

# Set RAPTOR_* variables with the content of ${CMAKE_FIND_PACKAGE_NAME}_(FOUND|...|VERSION)
# This needs to be done, because `find_package(RAPTOR)` might be called in any case-sensitive way and we want to
# guarantee that RAPTOR_* are always set.
foreach (package_var
         FOUND
         DIR
         ROOT
         CONFIG
         VERSION
         VERSION_MAJOR
         VERSION_MINOR
         VERSION_PATCH
         VERSION_TWEAK
         VERSION_COUNT
)
    set (RAPTOR_${package_var} "${${CMAKE_FIND_PACKAGE_NAME}_${package_var}}")
endforeach ()

# propagate RAPTOR_INCLUDE_DIR into RAPTOR_INCLUDE_DIRS
set (RAPTOR_INCLUDE_DIRS ${RAPTOR_INCLUDE_DIR} ${RAPTOR_DEPENDENCY_INCLUDE_DIRS})

# ----------------------------------------------------------------------------
# Export targets
# ----------------------------------------------------------------------------

if (RAPTOR_FOUND AND NOT TARGET raptor::raptor)
    separate_arguments (RAPTOR_CXX_FLAGS_LIST UNIX_COMMAND "${RAPTOR_CXX_FLAGS}")

    add_library (raptor_raptor INTERFACE)
    target_compile_definitions (raptor_raptor INTERFACE ${RAPTOR_DEFINITIONS})
    target_compile_options (raptor_raptor INTERFACE ${RAPTOR_CXX_FLAGS_LIST})
    target_link_options (raptor_raptor INTERFACE ${RAPTOR_CXX_FLAGS_LIST})
    target_link_libraries (raptor_raptor INTERFACE "${RAPTOR_LIBRARIES}" seqan::hibf sharg::sharg seqan3::seqan3)
    # include raptor/include/ as -I, because raptor should never produce warnings.
    target_include_directories (raptor_raptor INTERFACE "${RAPTOR_INCLUDE_DIR}")
    # include everything except raptor/include/ as -isystem, i.e.
    # a system header which suppresses warnings of external libraries.
    get_target_property (chopper_include_dir chopper_shared INCLUDE_DIRECTORIES)
    target_include_directories (raptor_raptor SYSTEM INTERFACE "${chopper_include_dir}")
    find_path (seqan3_test_include_dir
               NAMES seqan3/test/tmp_directory.hpp
               HINTS "${seqan3_SOURCE_DIR}/test/include"
    )
    target_include_directories (raptor_raptor SYSTEM INTERFACE "${seqan3_test_include_dir}")
    add_library (raptor::raptor ALIAS raptor_raptor)
endif ()

set (CMAKE_REQUIRED_QUIET ${CMAKE_REQUIRED_QUIET_SAVE})

if (RAPTOR_FIND_DEBUG)
    message ("Result for ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt")
    message ("")
    message ("  CMAKE_BUILD_TYPE            ${CMAKE_BUILD_TYPE}")
    message ("  CMAKE_SOURCE_DIR            ${CMAKE_SOURCE_DIR}")
    message ("  CMAKE_INCLUDE_PATH          ${CMAKE_INCLUDE_PATH}")
    message ("  RAPTOR_INCLUDE_DIR          ${RAPTOR_INCLUDE_DIR}")
    message ("")
    message ("  ${CMAKE_FIND_PACKAGE_NAME}_FOUND                ${${CMAKE_FIND_PACKAGE_NAME}_FOUND}")
    message ("  RAPTOR_HAS_ZLIB             ${ZLIB_FOUND}")
    message ("  RAPTOR_HAS_BZIP2            ${BZIP2_FOUND}")
    message ("")
    message ("  RAPTOR_INCLUDE_DIRS         ${RAPTOR_INCLUDE_DIRS}")
    message ("  RAPTOR_LIBRARIES            ${RAPTOR_LIBRARIES}")
    message ("  RAPTOR_DEFINITIONS          ${RAPTOR_DEFINITIONS}")
    message ("  RAPTOR_CXX_FLAGS            ${RAPTOR_CXX_FLAGS}")
    message ("")
    message ("  RAPTOR_VERSION              ${RAPTOR_VERSION}")
    message ("  RAPTOR_VERSION_MAJOR        ${RAPTOR_VERSION_MAJOR}")
    message ("  RAPTOR_VERSION_MINOR        ${RAPTOR_VERSION_MINOR}")
    message ("  RAPTOR_VERSION_PATCH        ${RAPTOR_VERSION_PATCH}")
endif ()

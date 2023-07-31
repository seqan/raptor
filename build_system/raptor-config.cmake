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
# RAPTOR_NO_BZIP2, RAPTOR_NO_CEREAL and RAPTOR_NO_LEMON respectively.
#
# If you wish to require the presence of ZLIB or BZip2, just check for the module before
# finding RAPTOR, e.g. "find_package (ZLIB REQUIRED)".
# If you wish to require the presence of CEREAL, you may define RAPTOR_CEREAL.
# If you wish to require the presence of LEMON, you may define RAPTOR_LEMON.
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

cmake_minimum_required (VERSION 3.4...3.12)

# ----------------------------------------------------------------------------
# Set initial variables
# ----------------------------------------------------------------------------

# make output globally quiet if required by find_package, this effects cmake functions like `check_*`
set (CMAKE_REQUIRED_QUIET_SAVE ${CMAKE_REQUIRED_QUIET})
set (CMAKE_REQUIRED_QUIET ${${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY})

# ----------------------------------------------------------------------------
# Greeter
# ----------------------------------------------------------------------------

string (ASCII 27 Esc)
set (ColourBold "${Esc}[1m")
set (ColourReset "${Esc}[m")

if (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
    message (STATUS "${ColourBold}Finding RAPTOR and checking requirements:${ColourReset}")
endif ()

# ----------------------------------------------------------------------------
# Includes
# ----------------------------------------------------------------------------

include (CheckIncludeFileCXX)
include (CheckCXXSourceCompiles)
include (FindPackageHandleStandardArgs)
include (CheckCXXCompilerFlag)

# ----------------------------------------------------------------------------
# Pretty printing and error handling
# ----------------------------------------------------------------------------

macro (raptor_config_print text)
    if (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
        message (STATUS "  ${text}")
    endif ()
endmacro ()

macro (raptor_config_error text)
    if (${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED)
        message (FATAL_ERROR ${text})
    else ()
        if (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
            message (WARNING ${text})
        endif ()
        return ()
    endif ()
endmacro ()

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
# Detect if we are a clone of repository and if yes auto-add submodules
# ----------------------------------------------------------------------------

if (RAPTOR_CLONE_DIR)
    raptor_config_print ("Detected as running from a repository checkout…")
endif ()

if (RAPTOR_SUBMODULES_DIR)
    file (GLOB
          submodules
          ${RAPTOR_SUBMODULES_DIR}/*/include
          ${RAPTOR_SUBMODULES_DIR}/*/test/include
          ${RAPTOR_SUBMODULES_DIR}/submodules/*/include
          ${RAPTOR_SUBMODULES_DIR}/*/src/include
          ${RAPTOR_SUBMODULES_DIR}/simde/simde
    )
    foreach (submodule ${submodules})
        if (IS_DIRECTORY ${submodule})
            raptor_config_print ("  …adding submodule include:  ${submodule}")
            set (RAPTOR_DEPENDENCY_INCLUDE_DIRS ${submodule} ${RAPTOR_DEPENDENCY_INCLUDE_DIRS})
        endif ()
    endforeach ()
    set (SHARG_HINT_TDL "${RAPTOR_SUBMODULES_DIR}/tool_description_lib")
endif ()

# ----------------------------------------------------------------------------
# Options for CheckCXXSourceCompiles
# ----------------------------------------------------------------------------

# deactivate messages in check_*
set (CMAKE_REQUIRED_QUIET 1)
# use global variables in Check* calls
set (CMAKE_REQUIRED_INCLUDES ${CMAKE_INCLUDE_PATH} ${RAPTOR_INCLUDE_DIR} ${RAPTOR_DEPENDENCY_INCLUDE_DIRS})
set (CMAKE_REQUIRED_FLAGS ${CMAKE_CXX_FLAGS})

# ----------------------------------------------------------------------------
# Force-deactivate optional dependencies
# ----------------------------------------------------------------------------

# These two are "opt-in", because detected by CMake
# If you want to force-require these, just do find_package (zlib REQUIRED) before find_package (raptor)
option (RAPTOR_NO_ZLIB "Don't use ZLIB, even if present." OFF)
option (RAPTOR_NO_BZIP2 "Don't use BZip2, even if present." OFF)

# ----------------------------------------------------------------------------
# Require C++20
# ----------------------------------------------------------------------------

set (CMAKE_REQUIRED_FLAGS_SAVE ${CMAKE_REQUIRED_FLAGS})

set (CXXSTD_TEST_SOURCE
     "#if !defined (__cplusplus) || (__cplusplus < 201709L)
    #error NOCXX20
    #endif
    int main() {}"
)

check_cxx_source_compiles ("${CXXSTD_TEST_SOURCE}" CXX20_BUILTIN)

if (CXX20_BUILTIN)
    raptor_config_print ("C++ Standard-20 support:    builtin")
else ()
    set (CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS_SAVE} -std=c++20")

    check_cxx_source_compiles ("${CXXSTD_TEST_SOURCE}" CXX20_FLAG)

    if (CXX20_FLAG)
        raptor_config_print ("C++ Standard-20 support:    via -std=c++20")
    else ()
        raptor_config_error ("RAPTOR requires C++20, but your compiler does not support it.")
    endif ()

    set (RAPTOR_CXX_FLAGS "${RAPTOR_CXX_FLAGS} -std=c++20")
endif ()

# ----------------------------------------------------------------------------
# Required: Sharg (with TDL for CWL/CTD support)
# ----------------------------------------------------------------------------
option (INSTALL_TDL "Enable installation of TDL." OFF)
find_package (Sharg QUIET REQUIRED HINTS ${RAPTOR_SUBMODULES_DIR}/sharg-parser/build_system)

# ----------------------------------------------------------------------------
# Optional: OpenMP
# ----------------------------------------------------------------------------

check_cxx_compiler_flag ("-fopenmp" RAPTOR_HAS_OPENMP)
if (RAPTOR_HAS_OPENMP)
    set (RAPTOR_CXX_FLAGS "${RAPTOR_CXX_FLAGS} -fopenmp")
    raptor_config_print ("OpenMP Support:             via -fopenmp")
else ()
    raptor_config_print ("OpenMP Support:             not found")
endif ()

check_cxx_compiler_flag ("-fopenmp-simd" RAPTOR_HAS_OPENMP_SIMD)
if (RAPTOR_HAS_OPENMP_SIMD)
    set (RAPTOR_CXX_FLAGS "${RAPTOR_CXX_FLAGS} -fopenmp-simd -DSIMDE_ENABLE_OPENMP")
    raptor_config_print ("SIMD-OpenMP Support:        via -fopenmp-simd")
else ()
    raptor_config_print ("SIMD-OpenMP Support:        not found")
endif ()

check_cxx_compiler_flag ("-Wno-psabi" RAPTOR_SUPPRESS_GCC4_ABI)
if (RAPTOR_SUPPRESS_GCC4_ABI)
    set (RAPTOR_CXX_FLAGS "${RAPTOR_CXX_FLAGS} -Wno-psabi")
    raptor_config_print ("Suppressing GCC 4 warnings: via -Wno-psabi")
else ()
    raptor_config_print ("Suppressing GCC 4 warnings: not found")
endif ()

# ----------------------------------------------------------------------------
# Optimizations
# ----------------------------------------------------------------------------

if ("${CMAKE_BUILD_TYPE}" MATCHES "Debug" OR "${CMAKE_BUILD_TYPE}" MATCHES "Coverage")
    set (RAPTOR_IS_DEBUG TRUE)
else ()
    set (RAPTOR_IS_DEBUG FALSE)
endif ()

option (RAPTOR_NATIVE_BUILD "Optimize build for current architecture." ON)
if (RAPTOR_IS_DEBUG)
    raptor_config_print ("Optimize build:             disabled")
elseif (RAPTOR_NATIVE_BUILD)
    set (RAPTOR_CXX_FLAGS "${RAPTOR_CXX_FLAGS} -march=native")
    raptor_config_print ("Optimize build:             via -march=native")
else ()
    check_cxx_compiler_flag ("-mpopcnt" RAPTOR_HAS_POPCNT)
    if (RAPTOR_HAS_POPCNT)
        set (RAPTOR_CXX_FLAGS "${RAPTOR_CXX_FLAGS} -mpopcnt")
        raptor_config_print ("Optimize build:             via -mpopcnt")
    else ()
        raptor_config_print ("Optimize build:             disabled")
    endif ()
endif ()

if (RAPTOR_IS_DEBUG)
    raptor_config_print ("Link Time Optimization:     disabled")
else ()
    if ("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
        set (RAPTOR_LTO_FLAGS "-flto=auto")
    else ()
        set (RAPTOR_LTO_FLAGS "-flto=auto -fno-fat-lto-objects")
    endif ()

    set (LTO_CMAKE_SOURCE
         "cmake_minimum_required(VERSION ${CMAKE_VERSION})\nproject(lto-test LANGUAGES CXX)
          cmake_policy(SET CMP0069 NEW)\nadd_library(foo foo.cpp)\nadd_executable(boo main.cpp)
          target_link_libraries(boo PUBLIC foo)"
    )
    set (LTO_FOO_CPP "int foo(){return 0x42;}")
    set (LTO_MAIN_CPP "int foo();int main(){return foo();}")
    set (testdir "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/_CMakeLTOTest-CXX")
    set (bindir "${testdir}/bin")
    set (srcdir "${testdir}/src")
    file (MAKE_DIRECTORY "${bindir}")
    file (MAKE_DIRECTORY "${srcdir}")
    file (WRITE "${srcdir}/foo.cpp" "${LTO_FOO_CPP}")
    file (WRITE "${srcdir}/main.cpp" "${LTO_MAIN_CPP}")
    file (WRITE "${srcdir}/CMakeLists.txt" "${LTO_CMAKE_SOURCE}")
    try_compile (RAPTOR_HAS_LTO "${bindir}"
                 "${srcdir}" "lto-test"
                 CMAKE_FLAGS "-DCMAKE_VERBOSE_MAKEFILE=ON" "-DCMAKE_CXX_FLAGS:STRING=${RAPTOR_LTO_FLAGS}"
                 OUTPUT_VARIABLE output
    )

    if (RAPTOR_HAS_LTO)
        raptor_config_print ("Link Time Optimization:     enabled")
        set (RAPTOR_CXX_FLAGS "${RAPTOR_CXX_FLAGS} ${RAPTOR_LTO_FLAGS}")
    else ()
        raptor_config_print ("Link Time Optimization:     not available")
    endif ()
endif ()

option (RAPTOR_STRIP_BINARY "Enable binary-stripping. Not supported on macOS." ON)
if ("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
    raptor_config_print ("Binary-stripping:           not supported")
elseif (RAPTOR_IS_DEBUG OR NOT RAPTOR_STRIP_BINARY)
    raptor_config_print ("Binary-stripping:           disabled")
else ()
    set (RAPTOR_CXX_FLAGS "${RAPTOR_CXX_FLAGS} -s")
    raptor_config_print ("Binary-stripping:           via -s")
endif ()

# ----------------------------------------------------------------------------
# thread support (pthread, windows threads)
# ----------------------------------------------------------------------------

set (THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package (Threads QUIET)

if (Threads_FOUND)
    set (RAPTOR_LIBRARIES ${RAPTOR_LIBRARIES} Threads::Threads)
    if ("${CMAKE_THREAD_LIBS_INIT}" STREQUAL "")
        raptor_config_print ("Thread support:             builtin")
    else ()
        raptor_config_print ("Thread support:             via ${CMAKE_THREAD_LIBS_INIT}")
    endif ()
else ()
    raptor_config_print ("Thread support:             not found")
endif ()

# ----------------------------------------------------------------------------
# ZLIB dependency
# ----------------------------------------------------------------------------

if (NOT RAPTOR_NO_ZLIB)
    find_package (ZLIB QUIET)
endif ()

if (ZLIB_FOUND)
    set (RAPTOR_LIBRARIES ${RAPTOR_LIBRARIES} ${ZLIB_LIBRARIES})
    set (RAPTOR_DEPENDENCY_INCLUDE_DIRS ${RAPTOR_DEPENDENCY_INCLUDE_DIRS} ${ZLIB_INCLUDE_DIRS})
    set (RAPTOR_DEFINITIONS ${RAPTOR_DEFINITIONS} "-DSEQAN3_HAS_ZLIB=1")
    raptor_config_print ("Optional dependency:        ZLIB-${ZLIB_VERSION_STRING} found")
else ()
    raptor_config_print ("Optional dependency:        ZLIB not found")
endif ()

# ----------------------------------------------------------------------------
# BZip2 dependency
# ----------------------------------------------------------------------------

if (NOT RAPTOR_NO_BZIP2)
    find_package (BZip2 QUIET)
endif ()

if (NOT ZLIB_FOUND AND BZIP2_FOUND)
    # NOTE (marehr): iostream_bzip2 uses the type `uInt`, which is defined by
    # `zlib`. Therefore, `bzip2` will cause a ton of errors without `zlib`.
    message (AUTHOR_WARNING "Disabling BZip2 [which was successfully found], "
                            "because ZLIB was not found BZip2 depends on ZLIB."
    )
    unset (BZIP2_FOUND)
endif ()

if (BZIP2_FOUND)
    set (RAPTOR_LIBRARIES ${RAPTOR_LIBRARIES} ${BZIP2_LIBRARIES})
    set (RAPTOR_DEPENDENCY_INCLUDE_DIRS ${RAPTOR_DEPENDENCY_INCLUDE_DIRS} ${BZIP2_INCLUDE_DIRS})
    set (RAPTOR_DEFINITIONS ${RAPTOR_DEFINITIONS} "-DSEQAN3_HAS_BZIP2=1")
    raptor_config_print ("Optional dependency:        BZip2-${BZIP2_VERSION_STRING} found")
else ()
    raptor_config_print ("Optional dependency:        BZip2 not found")
endif ()

# ----------------------------------------------------------------------------
# System dependencies
# ----------------------------------------------------------------------------

# librt
if ((${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    OR (${CMAKE_SYSTEM_NAME} STREQUAL "kFreeBSD")
    OR (${CMAKE_SYSTEM_NAME} STREQUAL "GNU")
)
    set (RAPTOR_LIBRARIES ${RAPTOR_LIBRARIES} rt)
endif ()

# libexecinfo -- implicit
check_include_file_cxx (execinfo.h _RAPTOR_HAVE_EXECINFO)
mark_as_advanced (_RAPTOR_HAVE_EXECINFO)
if (_RAPTOR_HAVE_EXECINFO)
    raptor_config_print ("Optional dependency:        libexecinfo found")
    if ((${CMAKE_SYSTEM_NAME} STREQUAL "FreeBSD") OR (${CMAKE_SYSTEM_NAME} STREQUAL "OpenBSD"))
        set (RAPTOR_LIBRARIES ${RAPTOR_LIBRARIES} execinfo elf)
    endif ()
else ()
    raptor_config_print ("Optional dependency:        libexecinfo not found")
endif ()

# ----------------------------------------------------------------------------
# Perform compilability test of platform.hpp (tests some requirements)
# ----------------------------------------------------------------------------

# set (CXXSTD_TEST_SOURCE
#      "#include <raptor/platform.hpp>
#      int main() {}")

# # using try_compile instead of check_cxx_source_compiles to capture output in case of failure
# file (WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.cxx" "${CXXSTD_TEST_SOURCE}\n")

# try_compile (RAPTOR_PLATFORM_TEST
#              ${CMAKE_BINARY_DIR}
#              ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.cxx
#              CMAKE_FLAGS         "-DCOMPILE_DEFINITIONS:STRING=${CMAKE_CXX_FLAGS} ${RAPTOR_CXX_FLAGS}"
#                                  "-DINCLUDE_DIRECTORIES:STRING=${CMAKE_INCLUDE_PATH};${RAPTOR_INCLUDE_DIR};${RAPTOR_DEPENDENCY_INCLUDE_DIRS}"
#              COMPILE_DEFINITIONS ${RAPTOR_DEFINITIONS}
#              LINK_LIBRARIES      ${RAPTOR_LIBRARIES}
#              OUTPUT_VARIABLE     RAPTOR_PLATFORM_TEST_OUTPUT)

# if (RAPTOR_PLATFORM_TEST)
#     raptor_config_print ("RAPTOR platform.hpp build:  passed.")
# else ()
#     raptor_config_error ("RAPTOR platform.hpp build:  failed!\n\
#                         ${RAPTOR_PLATFORM_TEST_OUTPUT}")
# endif ()

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
    target_link_libraries (raptor_raptor INTERFACE "${RAPTOR_LIBRARIES}" sharg::sharg)
    # include raptor/include/ as -I, because raptor should never produce warnings.
    target_include_directories (raptor_raptor INTERFACE "${RAPTOR_INCLUDE_DIR}")
    # include everything except raptor/include/ as -isystem, i.e.
    # a system header which suppresses warnings of external libraries.
    target_include_directories (raptor_raptor SYSTEM INTERFACE "${RAPTOR_DEPENDENCY_INCLUDE_DIRS}")
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

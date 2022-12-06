cmake_minimum_required (VERSION 3.16)

set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set (CMAKE_INSTALL_BINDIR "bin") # For chopper

# CCache

include (../test/cmake/raptor_require_ccache.cmake)
raptor_require_ccache ()

# Dependency: Chopper
include (FetchContent)

FetchContent_Declare (chopper_project URL "${RAPTOR_SUBMODULES_DIR}/chopper")

set (SEQAN3_SUBMODULES_DIR "${RAPTOR_SUBMODULES_DIR}")
set (CHOPPER_SUBMODULES_DIR
     "${RAPTOR_SUBMODULES_DIR}"
     CACHE STRING "Directory containing submodules."
)
set (CHOPPER_USE_CCACHE
     OFF
     CACHE BOOL "Use ccache if available."
)
set (CHOPPER_BUILD_DOC
     OFF
     CACHE BOOL "Build chopper documentation."
)
set (CHOPPER_BUILD_TEST
     OFF
     CACHE BOOL "Build chopper tests."
)
set (CMAKE_DISABLE_FIND_PACKAGE_Doxygen TRUE) # For chopper

FetchContent_MakeAvailable (chopper_project)

# Dependency: xxHash
FetchContent_Declare (raptor_xxhash_project URL "${RAPTOR_SUBMODULES_DIR}/xxHash" SOURCE_SUBDIR "cmake_unofficial")

set (XXHASH_BUILD_ENABLE_INLINE_API ON)
set (XXHASH_BUILD_XXHSUM OFF)
set (XXHASH_BUNDLED_MODE OFF)
set (BUILD_SHARED_LIBS OFF)

FetchContent_MakeAvailable (raptor_xxhash_project)

# Shared interface
if (NOT TARGET raptor_interface)
    add_library ("raptor_interface" INTERFACE)
    target_link_libraries ("raptor_interface" INTERFACE raptor::raptor)
    target_compile_options ("raptor_interface" INTERFACE "-pedantic" "-Wall" "-Wextra")
endif ()

# Raptor library
if (NOT TARGET raptor_lib)
    add_library ("raptor_lib" INTERFACE)
    target_link_libraries ("raptor_lib" INTERFACE "raptor_argument_parsing" "raptor_build" "raptor_build_hibf"
                                                  "raptor_search" "raptor_threshold" "raptor_upgrade" "raptor_layout"
    )
endif ()

# Raptor executable
if (NOT TARGET raptor)
    add_executable ("raptor" raptor.cpp)
    target_link_libraries ("raptor" PUBLIC "raptor_lib")
endif ()

install (TARGETS raptor chopper DESTINATION "${CMAKE_INSTALL_BINDIR}")

add_subdirectory (argument_parsing)
add_subdirectory (build)
add_subdirectory (layout)
add_subdirectory (search)
add_subdirectory (threshold)
add_subdirectory (upgrade)

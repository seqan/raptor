cmake_minimum_required (VERSION 3.18)

find_path (RAPTOR_CLONE_DIR
           NAMES build_system/raptor-config.cmake
           HINTS "${CMAKE_CURRENT_LIST_DIR}/../../.."
)

list (APPEND CMAKE_MODULE_PATH "${RAPTOR_CLONE_DIR}/build_system")

include (raptor-config-version)

project (raptor_utility_common LANGUAGES CXX C)

find_package (Raptor REQUIRED HINTS ${RAPTOR_CLONE_DIR}/build_system)

set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

if (NOT TARGET utility_common)
    add_library (utility_common INTERFACE)
    target_link_libraries ("utility_common" INTERFACE raptor::raptor)
    target_include_directories ("utility_common" INTERFACE "../include")
    target_compile_options ("utility_common" INTERFACE "-pedantic" "-Wall" "-Wextra")
endif ()

if (NOT TARGET hibf_sizes)
    add_executable ("hibf_sizes" hibf_sizes.cpp)
    target_link_libraries ("hibf_sizes" "utility_common")
endif ()

add_subdirectory (Genome_Biology)

add_subdirectory (iScience)

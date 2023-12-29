cmake_minimum_required (VERSION 3.21)

project (raptor_utility_common LANGUAGES CXX)

if (TARGET utility_common)
    return ()
endif ()

set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

get_filename_component (RAPTOR_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../.." ABSOLUTE)
option (RAPTOR_USE_CCACHE_IN_PARENT_PROJECT "" ON)
add_subdirectory ("${RAPTOR_ROOT_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/raptor")

add_library (utility_common INTERFACE)
target_link_libraries ("utility_common" INTERFACE raptor::raptor)
target_include_directories ("utility_common" INTERFACE "../include")

add_executable ("hibf_sizes" hibf_sizes.cpp)
target_link_libraries ("hibf_sizes" "utility_common")

add_subdirectory (Genome_Biology)

add_subdirectory (iScience)

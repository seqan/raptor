cmake_minimum_required (VERSION 3.8)

# Shared interface
add_library ("raptor_interface" INTERFACE)
target_link_libraries ("raptor_interface" INTERFACE seqan3::seqan3)
target_include_directories ("raptor_interface" INTERFACE ../include)
target_include_directories ("raptor_interface" INTERFACE ../lib/lemon/include)
target_include_directories ("raptor_interface" INTERFACE ../lib/robin-hood-hashing/src/include)

# Raptor library
add_library ("raptor_lib" INTERFACE)
target_link_libraries ("raptor_lib" INTERFACE "raptor_argument_parsing")
target_link_libraries ("raptor_lib" INTERFACE "raptor_build")
target_link_libraries ("raptor_lib" INTERFACE "raptor_search")
target_link_libraries ("raptor_lib" INTERFACE "raptor_upgrade")

# Raptor executable
add_executable ("raptor" raptor.cpp)
target_link_libraries ("raptor" PUBLIC "raptor_lib")

add_subdirectory (argument_parsing)
add_subdirectory (build)
add_subdirectory (search)
add_subdirectory (upgrade)

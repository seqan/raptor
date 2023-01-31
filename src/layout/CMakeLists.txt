cmake_minimum_required (VERSION 3.15)

add_library ("raptor_layout" STATIC raptor_layout.cpp)

target_link_libraries ("raptor_layout" PUBLIC "raptor_interface" "chopper_layout_lib" xxHash::xxhash)

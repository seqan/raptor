cmake_minimum_required (VERSION 3.15)

if (NOT TARGET raptor_layout)
    add_library ("raptor_layout" STATIC raptor_layout.cpp)
    target_compile_options (chopper_layout_lib PRIVATE "-w")
    target_link_libraries ("raptor_layout" PUBLIC "raptor_interface" "chopper_layout_lib" xxHash::xxhash)
endif ()

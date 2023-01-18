cmake_minimum_required (VERSION 3.18)

if (NOT TARGET raptor_prepare)
    add_library ("raptor_prepare" STATIC compute_minimiser.cpp)

    target_link_libraries ("raptor_prepare" PUBLIC "raptor_interface")
endif ()

cmake_minimum_required (VERSION 3.18)

if (NOT TARGET raptor_search)
    add_library ("raptor_search" STATIC raptor_search.cpp search_hibf.cpp search_ibf.cpp search_multiple.cpp)

    target_link_libraries ("raptor_search" PUBLIC "raptor_interface")
endif ()

cmake_minimum_required (VERSION 3.11)

add_library ("raptor_search" STATIC
             raptor_search.cpp
             search_hibf.cpp
             search_ibf.cpp
             search_multiple.cpp
             search_socks.cpp
)

target_link_libraries ("raptor_search" PUBLIC "raptor_interface")

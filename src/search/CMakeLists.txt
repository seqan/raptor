cmake_minimum_required (VERSION 3.21)

if (TARGET raptor_search)
    return ()
endif ()

add_library ("raptor_search" STATIC raptor_search.cpp search_hibf.cpp search_ibf.cpp search_partitioned_ibf.cpp)

target_link_libraries ("raptor_search" PUBLIC "raptor::raptor")

cmake_minimum_required (VERSION 3.8)

add_library ("raptor_precompute_threshold" STATIC precompute_threshold.cpp)
target_link_libraries ("raptor_precompute_threshold" PUBLIC "raptor_interface")

add_library ("raptor_search" STATIC raptor_search.cpp)
target_link_libraries ("raptor_search" PUBLIC "raptor_precompute_threshold" "raptor_search_detail")

add_subdirectory (detail)

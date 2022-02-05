cmake_minimum_required (VERSION 3.11)

add_library ("raptor_search" STATIC
             detail/one_indirect_error_model.cpp
             detail/multiple_error_model.cpp
             detail/pascal_row.cpp
             detail/one_error_model.cpp
             detail/precompute_correction.cpp
             detail/precompute_threshold.cpp
             raptor_search.cpp
             search_hibf.cpp
             search_ibf.cpp
             search_multiple.cpp
             search_socks.cpp
)

target_link_libraries ("raptor_search" PUBLIC "raptor_interface")

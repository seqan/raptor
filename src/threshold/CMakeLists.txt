cmake_minimum_required (VERSION 3.11)

add_library ("raptor_threshold" STATIC
             one_indirect_error_model.cpp
             multiple_error_model.cpp
             pascal_row.cpp
             one_error_model.cpp
             precompute_correction.cpp
             precompute_threshold.cpp
)

target_link_libraries ("raptor_threshold" PUBLIC "raptor_interface")

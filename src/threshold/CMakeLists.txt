cmake_minimum_required (VERSION 3.21)

if (TARGET raptor_threshold)
    return ()
endif ()

add_library ("raptor_threshold" STATIC
             multiple_error_model.cpp
             one_error_model.cpp
             one_indirect_error_model.cpp
             pascal_row.cpp
             precompute_correction.cpp
             precompute_threshold.cpp
             threshold.cpp
)

target_link_libraries ("raptor_threshold" PUBLIC "raptor::raptor")

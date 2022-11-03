cmake_minimum_required (VERSION 3.16)

if (NOT TARGET raptor_threshold)
    add_library ("raptor_threshold" STATIC multiple_error_model.cpp one_error_model.cpp one_indirect_error_model.cpp
                                           pascal_row.cpp precompute_correction.cpp precompute_threshold.cpp
    )

    target_link_libraries ("raptor_threshold" PUBLIC "raptor_interface")
endif ()

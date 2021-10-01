cmake_minimum_required (VERSION 3.8)

add_library ("raptor_search" STATIC detail/destroyed_indirectly_by_error.cpp
                                    detail/enumerate_all_errors.cpp
                                    detail/pascal_row.cpp
                                    detail/simple_model.cpp
                                    precompute_threshold.cpp
                                    raptor_search.cpp)

target_link_libraries ("raptor_search" PUBLIC "raptor_interface")

cmake_minimum_required (VERSION 3.8)

add_library ("raptor_search_pascal_row" STATIC pascal_row.cpp)
target_link_libraries ("raptor_search_pascal_row" PUBLIC "raptor_interface")

add_library ("raptor_search_simple_model" STATIC simple_model.cpp)
target_link_libraries ("raptor_search_simple_model" PUBLIC "raptor_search_pascal_row")

add_library ("raptor_search_enumerate_all_errors" STATIC enumerate_all_errors.cpp)
target_link_libraries ("raptor_search_enumerate_all_errors" PUBLIC "raptor_interface")

add_library ("raptor_search_destroyed_indirectly_by_error" STATIC destroyed_indirectly_by_error.cpp)
target_link_libraries ("raptor_search_destroyed_indirectly_by_error" PUBLIC "raptor_interface")

add_library ("raptor_search_detail" INTERFACE)
target_link_libraries ("raptor_search_detail" INTERFACE "raptor_search_pascal_row")
target_link_libraries ("raptor_search_detail" INTERFACE "raptor_search_simple_model")
target_link_libraries ("raptor_search_detail" INTERFACE "raptor_search_enumerate_all_errors")
target_link_libraries ("raptor_search_detail" INTERFACE "raptor_search_destroyed_indirectly_by_error")

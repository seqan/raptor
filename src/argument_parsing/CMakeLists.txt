cmake_minimum_required (VERSION 3.8)

add_library ("raptor_init_shared_meta" STATIC init_shared_meta.cpp)
target_compile_definitions ("raptor_init_shared_meta" PUBLIC "-DRAPTOR_VERSION=\"${CMAKE_PROJECT_VERSION} (${RAPTOR_COMMIT_HASH})\"")
target_compile_definitions ("raptor_init_shared_meta" PUBLIC "-DRAPTOR_DATE=\"${RAPTOR_COMMIT_DATE}\"")
target_link_libraries ("raptor_init_shared_meta" PUBLIC "raptor_interface")

add_library ("raptor_parse_bin_path" STATIC parse_bin_path.cpp)
target_link_libraries ("raptor_parse_bin_path" PUBLIC "raptor_interface")

add_library ("raptor_build_parsing" STATIC build_parsing.cpp)
target_link_libraries ("raptor_build_parsing" PUBLIC "raptor_init_shared_meta" "raptor_parse_bin_path")

add_library ("raptor_search_parsing" STATIC search_parsing.cpp)
target_link_libraries ("raptor_search_parsing" PUBLIC "raptor_init_shared_meta")

add_library ("raptor_upgrade_parsing" STATIC upgrade_parsing.cpp)
target_link_libraries ("raptor_upgrade_parsing" PUBLIC "raptor_init_shared_meta" "raptor_parse_bin_path")

add_library ("raptor_argument_parsing" INTERFACE)
target_link_libraries ("raptor_argument_parsing" INTERFACE "raptor_build_parsing")
target_link_libraries ("raptor_argument_parsing" INTERFACE "raptor_search_parsing")
target_link_libraries ("raptor_argument_parsing" INTERFACE "raptor_upgrade_parsing")

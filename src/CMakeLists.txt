cmake_minimum_required (VERSION 3.8)

# Shared interface
add_library ("${PROJECT_NAME}_interface" INTERFACE)
target_link_libraries ("${PROJECT_NAME}_interface" INTERFACE seqan3::seqan3)
target_include_directories ("${PROJECT_NAME}_interface" INTERFACE ../include)
target_include_directories ("${PROJECT_NAME}_interface" INTERFACE ../lib/robin-hood-hashing/src/include)

# Raptor build
add_library ("${PROJECT_NAME}_compute_minimiser_lib" STATIC build/compute_minimiser.cpp)
target_link_libraries ("${PROJECT_NAME}_compute_minimiser_lib" PUBLIC "${PROJECT_NAME}_interface")

add_library ("${PROJECT_NAME}_build_lib" STATIC build/raptor_build.cpp)
target_link_libraries ("${PROJECT_NAME}_build_lib" PUBLIC "${PROJECT_NAME}_compute_minimiser_lib")

# Raptor search
add_library ("${PROJECT_NAME}_search_helper_lib" STATIC search/detail/helper.cpp)
target_link_libraries ("${PROJECT_NAME}_search_helper_lib" PUBLIC "${PROJECT_NAME}_interface")

add_library ("${PROJECT_NAME}_minimiser_model_lib" STATIC search/minimiser_model.cpp)
target_link_libraries ("${PROJECT_NAME}_minimiser_model_lib" PUBLIC "${PROJECT_NAME}_search_helper_lib")

add_library ("${PROJECT_NAME}_simple_model_lib" STATIC search/compute_simple_model.cpp)
target_link_libraries ("${PROJECT_NAME}_simple_model_lib" PUBLIC "${PROJECT_NAME}_minimiser_model_lib")

add_library ("${PROJECT_NAME}_search_lib" STATIC search/raptor_search.cpp)
target_link_libraries ("${PROJECT_NAME}_search_lib" PUBLIC "${PROJECT_NAME}_simple_model_lib")

# Raptor upgrade
add_library ("${PROJECT_NAME}_upgrade_lib" STATIC upgrade/raptor_upgrade.cpp)
target_link_libraries ("${PROJECT_NAME}_upgrade_lib" PUBLIC "${PROJECT_NAME}_interface")

# Raptor argument parsing
add_library ("${PROJECT_NAME}_argument_parsing_init_shared_meta_lib" STATIC argument_parsing/init_shared_meta.cpp)
target_compile_definitions ("${PROJECT_NAME}_argument_parsing_init_shared_meta_lib" PUBLIC "-DRAPTOR_VERSION=\"${CMAKE_PROJECT_VERSION} (${RAPTOR_COMMIT_HASH})\"")
target_compile_definitions ("${PROJECT_NAME}_argument_parsing_init_shared_meta_lib" PUBLIC "-DRAPTOR_DATE=\"${RAPTOR_COMMIT_DATE}\"")
target_link_libraries ("${PROJECT_NAME}_argument_parsing_init_shared_meta_lib" PUBLIC "${PROJECT_NAME}_interface")

add_library ("${PROJECT_NAME}_argument_parsing_parse_bin_paths_lib" STATIC argument_parsing/parse_bin_paths.cpp)
target_link_libraries ("${PROJECT_NAME}_argument_parsing_parse_bin_paths_lib" PUBLIC "${PROJECT_NAME}_interface")

add_library ("${PROJECT_NAME}_argument_parsing_build_parsing_lib" STATIC argument_parsing/build_parsing.cpp)
target_link_libraries ("${PROJECT_NAME}_argument_parsing_build_parsing_lib" PUBLIC "${PROJECT_NAME}_argument_parsing_init_shared_meta_lib" "${PROJECT_NAME}_argument_parsing_parse_bin_paths_lib")

add_library ("${PROJECT_NAME}_argument_parsing_search_parsing_lib" STATIC argument_parsing/search_parsing.cpp)
target_link_libraries ("${PROJECT_NAME}_argument_parsing_search_parsing_lib" PUBLIC "${PROJECT_NAME}_argument_parsing_init_shared_meta_lib")

add_library ("${PROJECT_NAME}_argument_parsing_upgrade_parsing_lib" STATIC argument_parsing/upgrade_parsing.cpp)
target_link_libraries ("${PROJECT_NAME}_argument_parsing_upgrade_parsing_lib" PUBLIC "${PROJECT_NAME}_argument_parsing_init_shared_meta_lib" "${PROJECT_NAME}_argument_parsing_parse_bin_paths_lib")

# Raptor library
add_library ("${PROJECT_NAME}_lib" INTERFACE)
target_link_libraries ("${PROJECT_NAME}_lib" INTERFACE "${PROJECT_NAME}_argument_parsing_build_parsing_lib")
target_link_libraries ("${PROJECT_NAME}_lib" INTERFACE "${PROJECT_NAME}_argument_parsing_search_parsing_lib")
target_link_libraries ("${PROJECT_NAME}_lib" INTERFACE "${PROJECT_NAME}_argument_parsing_upgrade_parsing_lib")
target_link_libraries ("${PROJECT_NAME}_lib" INTERFACE "${PROJECT_NAME}_build_lib")
target_link_libraries ("${PROJECT_NAME}_lib" INTERFACE "${PROJECT_NAME}_search_lib")
target_link_libraries ("${PROJECT_NAME}_lib" INTERFACE "${PROJECT_NAME}_upgrade_lib")

# Raptor executable
add_executable ("${PROJECT_NAME}" raptor.cpp)
target_link_libraries ("${PROJECT_NAME}" PUBLIC "${PROJECT_NAME}_lib")

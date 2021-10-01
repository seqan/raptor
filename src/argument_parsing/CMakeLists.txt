cmake_minimum_required (VERSION 3.8)

add_library ("raptor_argument_parsing" STATIC build_parsing.cpp
                                              init_shared_meta.cpp
                                              init_shared_meta.cpp
                                              parse_bin_path.cpp
                                              search_parsing.cpp
                                              upgrade_parsing.cpp)

target_compile_definitions ("raptor_argument_parsing" PUBLIC
                            "-DRAPTOR_VERSION=\"${CMAKE_PROJECT_VERSION} (${RAPTOR_COMMIT_HASH})\""
                            "-DRAPTOR_DATE=\"${RAPTOR_COMMIT_DATE}\"")

target_link_libraries ("raptor_argument_parsing" PUBLIC "raptor_interface")


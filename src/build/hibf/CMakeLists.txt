cmake_minimum_required (VERSION 3.8)

add_library ("raptor_hibf_bin_size_in_bits" STATIC bin_size_in_bits.cpp)
target_link_libraries ("raptor_hibf_bin_size_in_bits" PUBLIC "raptor_interface")

add_library ("raptor_hibf_chopper_build" STATIC chopper_build.cpp)
target_link_libraries ("raptor_hibf_chopper_build" PUBLIC "raptor_hibf_create_ibfs_from_chopper_pack")

add_library ("raptor_hibf_compute_kmers" STATIC compute_kmers.cpp)
target_link_libraries ("raptor_hibf_compute_kmers" PUBLIC "raptor_interface")

add_library ("raptor_hibf_construct_ibf" STATIC construct_ibf.cpp)
target_link_libraries ("raptor_hibf_construct_ibf" PUBLIC "raptor_hibf_bin_size_in_bits")

add_library ("raptor_hibf_create_ibfs_from_chopper_pack" STATIC create_ibfs_from_chopper_pack.cpp)
target_link_libraries ("raptor_hibf_create_ibfs_from_chopper_pack" PUBLIC "raptor_hibf_read_chopper_pack_file" "raptor_hibf_hierarchical_build")

add_library ("raptor_hibf_hierarchical_build" STATIC hierarchical_build.cpp)
target_link_libraries ("raptor_hibf_hierarchical_build" PUBLIC "raptor_hibf_construct_ibf" "raptor_hibf_initialise_max_bin_kmers" "raptor_hibf_loop_over_children" "raptor_hibf_insert_into_ibf" "raptor_hibf_update_user_bins")

add_library ("raptor_hibf_initialise_max_bin_kmers" STATIC initialise_max_bin_kmers.cpp)
target_link_libraries ("raptor_hibf_initialise_max_bin_kmers" PUBLIC "raptor_hibf_compute_kmers")

add_library ("raptor_hibf_insert_into_ibf" STATIC insert_into_ibf.cpp)
target_link_libraries ("raptor_hibf_insert_into_ibf" PUBLIC "raptor_interface")

add_library ("raptor_hibf_loop_over_children" STATIC loop_over_children.cpp)
target_link_libraries ("raptor_hibf_loop_over_children" PUBLIC "raptor_interface")

add_library ("raptor_hibf_parse_chopper_pack_header" STATIC parse_chopper_pack_header.cpp)
target_link_libraries ("raptor_hibf_parse_chopper_pack_header" PUBLIC "raptor_interface")

add_library ("raptor_hibf_parse_chopper_pack_line" STATIC parse_chopper_pack_line.cpp)
target_link_libraries ("raptor_hibf_parse_chopper_pack_line" PUBLIC "raptor_interface")

add_library ("raptor_hibf_read_chopper_pack_file" STATIC read_chopper_pack_file.cpp)
target_link_libraries ("raptor_hibf_read_chopper_pack_file" PUBLIC "raptor_hibf_parse_chopper_pack_header" "raptor_hibf_parse_chopper_pack_line")

add_library ("raptor_hibf_update_user_bins" STATIC update_user_bins.cpp)
target_link_libraries ("raptor_hibf_update_user_bins" PUBLIC "raptor_interface")


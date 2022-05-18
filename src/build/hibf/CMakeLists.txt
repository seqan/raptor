cmake_minimum_required (VERSION 3.15)

add_library ("raptor_build_hibf" STATIC
             bin_size_in_bits.cpp
             chopper_build.cpp
             compute_kmers.cpp
             construct_ibf.cpp
             create_ibfs_from_chopper_pack.cpp
             hierarchical_build.cpp
             initialise_max_bin_kmers.cpp
             insert_into_ibf.cpp
             loop_over_children.cpp
             parse_chopper_pack_header.cpp
             parse_chopper_pack_line.cpp
             read_chopper_pack_file.cpp
             update_user_bins.cpp
)

target_link_libraries ("raptor_build_hibf" PUBLIC "raptor_interface")

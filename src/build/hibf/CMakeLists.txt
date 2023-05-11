cmake_minimum_required (VERSION 3.18)

if (NOT TARGET raptor_build_hibf)
    add_library ("raptor_build_hibf2" STATIC
                 bin_size_in_bits.cpp
                 chopper_build.cpp
                 compute_kmers.cpp
                 construct_ibf.cpp
                 create_ibfs_from_chopper_pack.cpp
                 hierarchical_build.cpp
                 initialise_build_tree.cpp
                 insert_into_ibf.cpp
                 loop_over_children.cpp
                 parse_chopper_pack_header.cpp
                 parse_chopper_pack_line.cpp
                 read_chopper_pack_file.cpp
                 update_content_node_data.cpp
                 update_header_node_data.cpp
                 update_parent_kmers.cpp
                 update_user_bins.cpp
    )
    target_link_libraries ("raptor_build_hibf2" PUBLIC "raptor_interface")

    # Version 1
    # add_library ("test_inserter" STATIC
    #              test_inserter.cpp
    # )
    # target_link_libraries ("test_inserter" PUBLIC "raptor_interface")

    # add_library ("raptor_build_hibf" INTERFACE)
    # target_link_libraries ("raptor_build_hibf" INTERFACE "raptor_build_hibf2" "test_inserter") # Order is important

    # Version 2
    add_library ("raptor_build_hibf" SHARED test_inserter.cpp)
    target_link_libraries ("raptor_build_hibf" PUBLIC "raptor_build_hibf2")
endif ()

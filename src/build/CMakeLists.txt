cmake_minimum_required (VERSION 3.18)

if (NOT TARGET raptor_build)
    add_library ("raptor_build" STATIC build_hibf.cpp build_ibf.cpp max_count_per_partition.cpp raptor_build.cpp)

    target_link_libraries ("raptor_build" PUBLIC "raptor_interface" "raptor_prepare" "seqan::hibf")

    # add_subdirectory (hibf)
endif ()

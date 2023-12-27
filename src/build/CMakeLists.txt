cmake_minimum_required (VERSION 3.21)

if (TARGET raptor_build)
    return ()
endif ()

add_library ("raptor_build" STATIC build_hibf.cpp build_ibf.cpp max_count_per_partition.cpp raptor_build.cpp)

target_link_libraries ("raptor_build" PUBLIC "raptor::raptor" "raptor_prepare" "seqan::hibf")

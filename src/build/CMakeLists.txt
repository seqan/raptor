cmake_minimum_required (VERSION 3.8)

add_library ("raptor_compute_minimiser" STATIC compute_minimiser.cpp)
target_link_libraries ("raptor_compute_minimiser" PUBLIC "raptor_interface")

add_library ("raptor_build_from_minimiser" STATIC build_from_minimiser.cpp)
target_link_libraries ("raptor_build_from_minimiser" PUBLIC "raptor_interface")

add_library ("raptor_build" STATIC raptor_build.cpp)
target_link_libraries ("raptor_build" PUBLIC "raptor_compute_minimiser" "raptor_build_from_minimiser" "raptor_hibf_chopper_build")

add_subdirectory (hibf)

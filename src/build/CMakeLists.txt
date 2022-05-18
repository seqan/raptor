cmake_minimum_required (VERSION 3.15)

add_library ("raptor_build" STATIC build_from_minimiser.cpp compute_minimiser.cpp raptor_build.cpp)

target_link_libraries ("raptor_build" PUBLIC "raptor_interface")

add_subdirectory (hibf)

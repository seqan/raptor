cmake_minimum_required (VERSION 3.25)

if (TARGET raptor_prepare)
    return ()
endif ()

add_library ("raptor_prepare" STATIC compute_minimiser.cpp)

target_link_libraries ("raptor_prepare" PUBLIC "raptor::raptor")

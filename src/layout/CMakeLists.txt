cmake_minimum_required (VERSION 3.15)

if (TARGET raptor_layout)
    return ()
endif ()

add_library ("raptor_layout" STATIC raptor_layout.cpp)
target_link_libraries ("raptor_layout" PUBLIC "raptor::raptor" "chopper_lib")

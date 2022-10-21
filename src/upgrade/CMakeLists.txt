cmake_minimum_required (VERSION 3.15)

if (NOT TARGET raptor_upgrade)
    add_library ("raptor_upgrade" STATIC raptor_upgrade.cpp)
    target_link_libraries ("raptor_upgrade" PUBLIC "raptor_interface")
endif ()

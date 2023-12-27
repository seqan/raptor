cmake_minimum_required (VERSION 3.21)

if (TARGET raptor_upgrade)
    return ()
endif ()

add_library ("raptor_upgrade" STATIC raptor_upgrade.cpp)
target_link_libraries ("raptor_upgrade" PUBLIC "raptor::raptor")

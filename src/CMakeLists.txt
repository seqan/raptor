cmake_minimum_required (VERSION 3.21)

if (TARGET raptor)
    return ()
endif ()

# Raptor library
add_library ("raptor_lib" INTERFACE)
target_link_libraries ("raptor_lib"
                       INTERFACE "raptor_argument_parsing"
                                 "raptor_build"
                                 "raptor_prepare"
                                 "raptor_search"
                                 "raptor_threshold"
                                 "raptor_upgrade"
                                 "raptor_layout"
)

# Raptor executable
add_executable ("raptor" raptor.cpp)
target_link_libraries ("raptor" PUBLIC "raptor_lib")

add_subdirectory (argument_parsing)
add_subdirectory (build)
add_subdirectory (layout)
add_subdirectory (search)
add_subdirectory (prepare)
add_subdirectory (threshold)
add_subdirectory (upgrade)

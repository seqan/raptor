cmake_minimum_required (VERSION 3.18)

set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set (CMAKE_INSTALL_BINDIR "bin") # For chopper

# Shared interface
if (NOT TARGET raptor_interface)
    add_library ("raptor_interface" INTERFACE)
    target_link_libraries ("raptor_interface" INTERFACE raptor::raptor)
    target_compile_options ("raptor_interface" INTERFACE "-pedantic" "-Wall" "-Wextra")
endif ()

# Raptor library
if (NOT TARGET raptor_lib)
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
endif ()

# Raptor executable
if (NOT TARGET raptor)
    add_executable ("raptor" raptor.cpp)
    target_link_libraries ("raptor" PUBLIC "raptor_lib")
endif ()

install (TARGETS raptor RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")

add_subdirectory (argument_parsing)
add_subdirectory (build)
add_subdirectory (layout)
add_subdirectory (search)
add_subdirectory (prepare)
add_subdirectory (threshold)
add_subdirectory (upgrade)

cmake_minimum_required (VERSION 3.15)

# Shared interface
add_library ("raptor_interface" INTERFACE)
target_link_libraries ("raptor_interface" INTERFACE seqan3::seqan3 chopper_lib)
target_include_directories ("raptor_interface" INTERFACE ../include)
target_include_directories ("raptor_interface" SYSTEM INTERFACE ${RAPTOR_SUBMODULES_DIR}/chopper/include)
target_include_directories ("raptor_interface" SYSTEM INTERFACE ${RAPTOR_SUBMODULES_DIR}/robin-hood-hashing/src/include)
target_compile_options ("raptor_interface" INTERFACE "-pedantic" "-Wall" "-Wextra")
add_dependencies (raptor_interface chopper_lib)

# Raptor library
add_library ("raptor_lib" INTERFACE)
target_link_libraries ("raptor_lib" INTERFACE "raptor_argument_parsing" "raptor_build" "raptor_build_hibf"
                                              "raptor_search" "raptor_threshold" "raptor_upgrade"
)

# Raptor executable
add_executable ("raptor" raptor.cpp)
target_link_libraries ("raptor" PUBLIC "raptor_lib")

install (TARGETS raptor DESTINATION bin)
if (NOT CMAKE_VERSION VERSION_LESS 3.21) # cmake >= 3.21
    install (IMPORTED_RUNTIME_ARTIFACTS chopper DESTINATION bin)
else ()
    install (PROGRAMS "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/chopper" DESTINATION bin)
endif ()

add_subdirectory (argument_parsing)
add_subdirectory (build)
add_subdirectory (search)
add_subdirectory (threshold)
add_subdirectory (upgrade)

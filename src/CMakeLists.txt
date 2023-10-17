cmake_minimum_required (VERSION 3.18)

set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set (CMAKE_INSTALL_BINDIR "bin") # For chopper

include (FetchContent)

# Dependency: HIBF lib
set (HIBF_SUBMODULES_DIR
     "${RAPTOR_SUBMODULES_DIR}"
     CACHE STRING "Directory containing submodules."
)

if (NOT TARGET seqan::hibf)
    FetchContent_Declare (raptor_hibf_project SOURCE_DIR "${RAPTOR_SUBMODULES_DIR}/hibf")
    option (INSTALL_HIBF "Enable installation of HIBF." OFF)
    option (HIBF_NATIVE_BUILD "Optimize build for current architecture." ${RAPTOR_NATIVE_BUILD})
    set (STORED_CMAKE_MESSAGE_LOG_LEVEL "${CMAKE_MESSAGE_LOG_LEVEL}")
    set (CMAKE_MESSAGE_LOG_LEVEL "ERROR")
    FetchContent_MakeAvailable (raptor_hibf_project)
    set (CMAKE_MESSAGE_LOG_LEVEL "${STORED_CMAKE_MESSAGE_LOG_LEVEL}")
    unset (STORED_CMAKE_MESSAGE_LOG_LEVEL)
endif ()

# Dependency: Chopper
FetchContent_Declare (raptor_chopper_project SOURCE_DIR "${RAPTOR_SUBMODULES_DIR}/chopper")

set (SEQAN3_SUBMODULES_DIR "${RAPTOR_SUBMODULES_DIR}")
set (CHOPPER_SUBMODULES_DIR
     "${RAPTOR_SUBMODULES_DIR}"
     CACHE STRING "Directory containing submodules."
)
set (CHOPPER_USE_CCACHE
     OFF
     CACHE BOOL "Use ccache if available."
)
set (CHOPPER_BUILD_DOC
     OFF
     CACHE BOOL "Build chopper documentation."
)
set (CHOPPER_BUILD_TEST
     OFF
     CACHE BOOL "Build chopper tests."
)
set (CHOPPER_INSTALL
     OFF
     CACHE BOOL "Install Chopper."
)

set (STORED_CMAKE_MESSAGE_LOG_LEVEL "${CMAKE_MESSAGE_LOG_LEVEL}")
set (CMAKE_MESSAGE_LOG_LEVEL "ERROR")
FetchContent_MakeAvailable (raptor_chopper_project)
set_target_properties (chopper PROPERTIES EXCLUDE_FROM_ALL TRUE)
set (CMAKE_MESSAGE_LOG_LEVEL "${STORED_CMAKE_MESSAGE_LOG_LEVEL}")
unset (STORED_CMAKE_MESSAGE_LOG_LEVEL)

# Shared interface
if (NOT TARGET raptor_interface)
    add_library ("raptor_interface" INTERFACE)
    target_link_libraries ("raptor_interface" INTERFACE raptor::raptor)
    target_compile_options ("raptor_interface" INTERFACE "-pedantic" "-Wall" "-Wextra")
    target_compile_definitions ("raptor_interface" INTERFACE "RAPTOR_OLD_HIBF")
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

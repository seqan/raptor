cmake_minimum_required (VERSION 3.8)

# Register how to download Googletest.
include (ExternalProject)
ExternalProject_Add (googletest
                     GIT_REPOSITORY    "https://github.com/google/googletest.git"
                     GIT_TAG           "master"
                     GIT_SHALLOW
                     PREFIX            "${CMAKE_CURRENT_BINARY_DIR}/googletest"
                     UPDATE_COMMAND    ""   # omit update step
                     INSTALL_COMMAND   "")  # do not install

# The 4 library targets of googletest are combined in the gtest_all interface.
add_library (gtest_all INTERFACE)
ExternalProject_Get_Property (googletest binary_dir)
foreach (target "gtest" "gmock" "gtest_main" "gmock_main")
    # Import the library from the googletest build directory.
    add_library (${target} UNKNOWN IMPORTED)
    set_target_properties (${target} PROPERTIES
                           IMPORTED_LOCATION "${binary_dir}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}${target}.a")
    # Require googletest to be downloaded before the library is created.
    add_dependencies (${target} googletest)
    # Add the library to the 'gtest_all' interface target.
    target_link_libraries (gtest_all INTERFACE ${target})
endforeach ()

# Add the include directories to the 'gtest_all' interface target.
ExternalProject_Get_Property (googletest source_dir)
target_include_directories (gtest_all INTERFACE "${source_dir}/googletest/include")
target_include_directories (gtest_all INTERFACE "${source_dir}/googlemock/include")

cmake_minimum_required (VERSION 3.8)

# Set the project specific configuration settings to configure and build the gtest target.
set (gtest_project_args ${APP_TEMPLATE_EXTERNAL_PROJECT_CMAKE_ARGS})

# Force that libraries are installed to `lib/`, because GNUInstallDirs might install it into `lib64/`.
list (APPEND gtest_project_args "-DCMAKE_INSTALL_LIBDIR=${PROJECT_BINARY_DIR}/lib/")

# Register how to download Googletest.
include (ExternalProject)
ExternalProject_Add (googletest
                     GIT_REPOSITORY    "https://github.com/google/googletest.git"
                     GIT_TAG           "master"
                     GIT_SHALLOW
                     PREFIX            "${CMAKE_CURRENT_BINARY_DIR}/googletest"
                     CMAKE_ARGS        "${gtest_project_args}"
                     UPDATE_COMMAND    ""   # omit update step
                     INSTALL_COMMAND   "")  # do not install

unset (gtest_project_args)

# The 4 library targets of googletest are combined in the gtest_all interface.
add_library (gtest_all INTERFACE)
ExternalProject_Get_Property (googletest binary_dir)
foreach (target "gtest" "gmock" "gtest_main" "gmock_main")
    # Import the library from the googletest build directory.
    add_library (${target} UNKNOWN IMPORTED)

    # Google test adds 'd' to its targets in debug build mode.
    # This makes sure, that the target is linked against the correct compatible library.
    set (target_suffix "")
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        set (target_suffix "d")
    endif ()

    # Define the proper target location.
    set (target_location
        "${binary_dir}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}${target}${target_suffix}${CMAKE_STATIC_LIBRARY_SUFFIX}")

    set_target_properties (${target} PROPERTIES IMPORTED_LOCATION "${target_location}")
    # Require googletest to be downloaded before the library is created.
    add_dependencies (${target} googletest)
    # Add the library to the 'gtest_all' interface target.
    target_link_libraries (gtest_all INTERFACE ${target})

    unset (target_suffix)
    unset (target_location)
endforeach ()

# Add the include directories to the 'gtest_all' interface target.
ExternalProject_Get_Property (googletest source_dir)
target_include_directories (gtest_all INTERFACE "${source_dir}/googletest/include")
target_include_directories (gtest_all INTERFACE "${source_dir}/googlemock/include")

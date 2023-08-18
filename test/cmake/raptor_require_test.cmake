# --------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
# --------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.18)

# Exposes the google-test targets `gtest` and `gtest_main`.
# CMake 3.24: https://cmake.org/cmake/help/latest/module/FetchContent.html#variable:FETCHCONTENT_TRY_FIND_PACKAGE_MODE
macro (raptor_require_test)
    enable_testing ()

    set (RAPTOR_GTEST_TAG "v1.13.0")

    find_package (GTest QUIET)

    # Also ensure that Google Test if fetched for the latest library cron, which sets the tag to "main".
    if (NOT GTest_FOUND OR "${RAPTOR_GTEST_TAG}" STREQUAL "main")
        message (STATUS "  Fetching Google Test ${RAPTOR_GTEST_TAG}")

        include (FetchContent)
        FetchContent_Declare (gtest_fetch_content
                              GIT_REPOSITORY "https://github.com/google/googletest.git"
                              GIT_TAG "${RAPTOR_GTEST_TAG}"
        )
        option (BUILD_GMOCK "" OFF)
        option (INSTALL_GTEST "" OFF)
        set (STORED_CMAKE_MESSAGE_LOG_LEVEL "${CMAKE_MESSAGE_LOG_LEVEL}")
        set (CMAKE_MESSAGE_LOG_LEVEL "ERROR")
        FetchContent_MakeAvailable (gtest_fetch_content)
        set (CMAKE_MESSAGE_LOG_LEVEL "${STORED_CMAKE_MESSAGE_LOG_LEVEL}")
        unset (STORED_CMAKE_MESSAGE_LOG_LEVEL)
    else ()
        message (STATUS "  Test dependency:            Google Test ${GTest_VERSION} found.")
    endif ()

    if (NOT TARGET gtest_build)
        add_custom_target (gtest_build DEPENDS gtest_main gtest)
    endif ()

endmacro ()

# --------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
# --------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.18)

# Exposes the google-benchmark target `benchmark` and `benchmark_main`.
# CMake 3.24: https://cmake.org/cmake/help/latest/module/FetchContent.html#variable:FETCHCONTENT_TRY_FIND_PACKAGE_MODE
macro (raptor_require_benchmark)
    enable_testing ()

    set (RAPTOR_BENCHMARK_TAG "v1.8.2")

    find_package (benchmark QUIET)

    # Also ensure that Google Benchmark if fetched for the latest library cron, which sets the tag to "main".
    if (NOT benchmark_FOUND OR "${RAPTOR_BENCHMARK_TAG}" STREQUAL "main")
        message (STATUS "  Fetching Google Benchmark ${RAPTOR_BENCHMARK_TAG}")

        include (FetchContent)
        FetchContent_Declare (gbenchmark_fetch_content
                              GIT_REPOSITORY "https://github.com/google/benchmark.git"
                              GIT_TAG "${RAPTOR_BENCHMARK_TAG}"
        )
        option (BENCHMARK_ENABLE_TESTING "" OFF)
        option (BENCHMARK_ENABLE_WERROR "" OFF) # Does not apply to Debug builds.
        option (BENCHMARK_ENABLE_INSTALL "" OFF)
        set (STORED_CMAKE_MESSAGE_LOG_LEVEL "${CMAKE_MESSAGE_LOG_LEVEL}")
        set (CMAKE_MESSAGE_LOG_LEVEL "ERROR")
        FetchContent_MakeAvailable (gbenchmark_fetch_content)
        set (CMAKE_MESSAGE_LOG_LEVEL "${STORED_CMAKE_MESSAGE_LOG_LEVEL}")
        unset (STORED_CMAKE_MESSAGE_LOG_LEVEL)
    else ()
        message (STATUS "  Test dependency:            Google Benchmark ${benchmark_VERSION} found.")
    endif ()

    # NOTE: google benchmark's CMakeLists.txt already defines Shlwapi
    if (NOT TARGET gbenchmark_build)
        add_custom_target (gbenchmark_build DEPENDS benchmark_main benchmark)
    endif ()

endmacro ()

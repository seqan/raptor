# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.15)

# Exposes the google-benchmark target `gbenchmark`.
macro (raptor_require_benchmark)
    enable_testing ()

    set (RAPTOR_GBENCHMARK_GIT_TAG "37be1e8252527229cccad9f097afe68572f3c08a" CACHE STRING "googlebenchmark commit to use")

    message (STATUS "Fetch Google Benchmark:")

    include (FetchContent)
    FetchContent_Declare (
        gbenchmark_fetch_content
        GIT_REPOSITORY "https://github.com/google/benchmark.git"
        GIT_TAG "${RAPTOR_GBENCHMARK_GIT_TAG}")
    option (BENCHMARK_ENABLE_TESTING "" OFF)
    option (BENCHMARK_ENABLE_WERROR "" OFF) # Does not apply to Debug builds.
    FetchContent_MakeAvailable (gbenchmark_fetch_content)

    # NOTE: google benchmark's CMakeLists.txt already defines Shlwapi
    add_library (gbenchmark ALIAS benchmark_main)

    add_custom_target (gbenchmark_build DEPENDS gbenchmark)
    if (NOT TARGET gtest_build)
        add_custom_target (gbenchmark_build DEPENDS gbenchmark)
        target_compile_options ("gbenchmark" PUBLIC "-w")
    endif ()
endmacro ()

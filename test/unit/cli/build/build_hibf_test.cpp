// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h> // for Message, TestParamInfo, TestPartResult, AssertionResult

#include <algorithm>  // for max
#include <cstddef>    // for size_t
#include <filesystem> // for path, exists
#include <ranges>     // for operator|, operator==, views
#include <string>     // for basic_string, operator+, string, to_string, allocator
#include <tuple>      // for get, tuple

#include <seqan3/test/pretty_printing.hpp>             // for PrintTo
#include <seqan3/utility/container/dynamic_bitset.hpp> // for operator==

#include <raptor/index.hpp>         // for hibf
#include <raptor/test/cli_test.hpp> // for RAPTOR_ASSERT_ZERO_EXIT, raptor_base, raptor_base::st...

struct build_hibf : public raptor_base, public testing::WithParamInterface<std::tuple<size_t, size_t, bool>>
{};

TEST_P(build_hibf, with_file)
{
    auto const [number_of_repeated_bins, window_size, run_parallel_tmp] = GetParam();
    bool const run_parallel = run_parallel_tmp && number_of_repeated_bins >= 32;

    cli_test_result const result = execute_app("raptor",
                                               "build",
                                               "--kmer 19",
                                               "--window",
                                               std::to_string(window_size),
                                               "--threads",
                                               run_parallel ? "2" : "1",
                                               "--output raptor.index",
                                               "--quiet",
                                               "--input",
                                               layout_path(number_of_repeated_bins));
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);

    compare_index<raptor::index_structure::hibf>(ibf_path(number_of_repeated_bins, window_size, is_hibf::yes),
                                                 "raptor.index");
}

TEST_P(build_hibf, with_shape)
{
    auto const [number_of_repeated_bins, window_size, run_parallel_tmp] = GetParam();
    bool const run_parallel = run_parallel_tmp && number_of_repeated_bins >= 32;

    cli_test_result const result = execute_app("raptor",
                                               "build",
                                               "--shape 1111111111111111111",
                                               "--window",
                                               std::to_string(window_size),
                                               "--threads",
                                               run_parallel ? "2" : "1",
                                               "--output raptor.index",
                                               "--quiet",
                                               "--input",
                                               layout_path(number_of_repeated_bins));
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);

    compare_index<raptor::index_structure::hibf>(ibf_path(number_of_repeated_bins, window_size, is_hibf::yes),
                                                 "raptor.index");
}

INSTANTIATE_TEST_SUITE_P(build_hibf_suite,
                         build_hibf,
                         testing::Combine(testing::Values(0, 16, 32),
                                          testing::Values(19, 23),
                                          testing::Values(true, false)),
                         [](testing::TestParamInfo<build_hibf::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::max<int>(1, std::get<0>(info.param) * 4)) + "_bins_"
                                              + std::to_string(std::get<1>(info.param)) + "_window_"
                                              + (std::get<2>(info.param) ? "parallel" : "serial");
                             return name;
                         });

TEST_F(build_hibf, three_levels)
{
    cli_test_result const result = execute_app("raptor",
                                               "build",
                                               "--kmer 19",
                                               "--window 19",
                                               "--threads 1",
                                               "--output raptor.index",
                                               "--quiet",
                                               "--input",
                                               data("three_levels.layout"));
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);

    compare_index<raptor::index_structure::hibf>(data("three_levels.hibf"), "raptor.index");
}

TEST_F(build_hibf, verbose)
{
    cli_test_result const result = execute_app("raptor",
                                               "build",
                                               "--kmer 19",
                                               "--window 19",
                                               "--threads 1",
                                               "--output raptor.index",
                                               "--timing-output raptor.time",
                                               "--input",
                                               data("three_levels.layout"));
    EXPECT_EQ(result.out, std::string{});
    EXPECT_TRUE(result.err.starts_with("============= Timings ============="));
    EXPECT_TRUE(std::filesystem::exists("raptor.time"));
    RAPTOR_ASSERT_ZERO_EXIT(result);

    compare_index<raptor::index_structure::hibf>(data("three_levels.hibf"), "raptor.index");
}

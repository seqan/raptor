// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h> // for Message, TestParamInfo, AssertionResult, Values, Test...

#include <algorithm>  // for max
#include <cstddef>    // for size_t
#include <filesystem> // for path, exists
#include <fstream>    // for basic_ofstream, operator<<, basic_ostream, ofstream
#include <ranges>     // for join_view
#include <string>     // for basic_string, operator+, char_traits, to_string, string
#include <tuple>      // for get, tuple

#include <seqan3/test/pretty_printing.hpp>             // for PrintTo
#include <seqan3/utility/container/dynamic_bitset.hpp> // for operator==

#include <hibf/contrib/std/pair.hpp>         // for get
#include <hibf/interleaved_bloom_filter.hpp> // for interleaved_bloom_filter

#include <raptor/test/cli_test.hpp> // for RAPTOR_ASSERT_ZERO_EXIT, raptor_base

struct build_ibf : public raptor_base, public testing::WithParamInterface<std::tuple<size_t, size_t, bool>>
{};

TEST_P(build_ibf, with_file)
{
    auto const [number_of_repeated_bins, window_size, run_parallel_tmp] = GetParam();
    bool const run_parallel = run_parallel_tmp && number_of_repeated_bins >= 32;

    { // generate input file
        std::ofstream file{"raptor_cli_test.txt"};
        for (auto && file_path : get_repeated_bins(number_of_repeated_bins))
            file << file_path << '\n';
        file << '\n';
    }

    cli_test_result const result = execute_app("raptor",
                                               "build",
                                               "--kmer 19",
                                               "--window ",
                                               std::to_string(window_size),
                                               "--threads ",
                                               run_parallel ? "2" : "1",
                                               "--output raptor.index",
                                               "--quiet",
                                               "--input",
                                               "raptor_cli_test.txt");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);

    compare_index(ibf_path(number_of_repeated_bins, window_size), "raptor.index");
}

TEST_P(build_ibf, with_shape)
{
    auto const [number_of_repeated_bins, window_size, run_parallel_tmp] = GetParam();
    bool const run_parallel = run_parallel_tmp && number_of_repeated_bins >= 32;

    { // generate input file
        std::ofstream file{"raptor_cli_test.txt"};
        for (auto && file_path : get_repeated_bins(number_of_repeated_bins))
            file << file_path << '\n';
        file << '\n';
    }

    cli_test_result const result = execute_app("raptor",
                                               "build",
                                               "--shape 1111111111111111111",
                                               "--window ",
                                               std::to_string(window_size),
                                               "--threads ",
                                               run_parallel ? "2" : "1",
                                               "--output raptor.index",
                                               "--quiet",
                                               "--input",
                                               "raptor_cli_test.txt");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);

    compare_index(ibf_path(number_of_repeated_bins, window_size), "raptor.index");
}

INSTANTIATE_TEST_SUITE_P(build_ibf_suite,
                         build_ibf,
                         testing::Combine(testing::Values(0, 16, 32),
                                          testing::Values(19, 23),
                                          testing::Values(true, false)),
                         [](testing::TestParamInfo<build_ibf::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::max<int>(1, std::get<0>(info.param) * 4)) + "_bins_"
                                              + std::to_string(std::get<1>(info.param)) + "_window_"
                                              + (std::get<2>(info.param) ? "parallel" : "serial");
                             return name;
                         });

TEST_F(build_ibf, verbose)
{
    { // generate input file
        std::ofstream file{"raptor_cli_test.txt"};
        for (auto && file_path : get_repeated_bins(16u))
            file << file_path << '\n';
        file << '\n';
    }

    cli_test_result const result = execute_app("raptor",
                                               "build",
                                               "--kmer 19",
                                               "--window 19",
                                               "--threads 1",
                                               "--output raptor.index",
                                               "--timing-output raptor.time",
                                               "--input",
                                               "raptor_cli_test.txt");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_TRUE(result.err.starts_with("============= Timings ============="));
    EXPECT_TRUE(std::filesystem::exists("raptor.time"));
    RAPTOR_ASSERT_ZERO_EXIT(result);

    compare_index(ibf_path(16, 19), "raptor.index");
}

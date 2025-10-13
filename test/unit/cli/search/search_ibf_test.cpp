// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h> // for Message, TestParamInfo, TestPartResult, CmpHelperEQ, CmpHelperEQF...

#include <algorithm>   // for max
#include <cstddef>     // for size_t
#include <cstdint>     // for uint32_t, uint8_t
#include <filesystem>  // for exists, path
#include <string>      // for basic_string, operator+, to_string, string, allocator, char_traits
#include <string_view> // for basic_string_view
#include <tuple>       // for get, tuple

#include <seqan3/test/pretty_printing.hpp> // for PrintTo

#include <raptor/test/cli_test.hpp> // for RAPTOR_ASSERT_ZERO_EXIT, raptor_base, raptor_base::strong_bool::yes

struct search_ibf : public raptor_base, public testing::WithParamInterface<std::tuple<size_t, size_t, size_t>>
{};

TEST_P(search_ibf, with_error)
{
    auto const [number_of_repeated_bins, window_size, number_of_errors] = GetParam();

    cli_test_result const result = execute_app("raptor",
                                               "search",
                                               "--output search.out",
                                               "--error ",
                                               std::to_string(number_of_errors),
                                               "--p_max 0.4",
                                               "--index ",
                                               ibf_path(number_of_repeated_bins, window_size),
                                               "--quiet",
                                               "--query ",
                                               data("query.fq"));
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);

    compare_search(number_of_repeated_bins, number_of_errors, "search.out");
}

TEST_P(search_ibf, threshold)
{
    auto const [number_of_repeated_bins, window_size, number_of_errors] = GetParam();

    cli_test_result const result = execute_app("raptor",
                                               "search",
                                               "--output search.out",
                                               "--threshold 0.50",
                                               "--index ",
                                               ibf_path(number_of_repeated_bins, window_size),
                                               "--quiet",
                                               "--query ",
                                               data("query.fq"));
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);

    compare_search(number_of_repeated_bins, 1 /* Always finds everything */, "search.out");
}

TEST_P(search_ibf, no_hits)
{
    auto const [number_of_repeated_bins, window_size, number_of_errors] = GetParam();

    cli_test_result const result = execute_app("raptor",
                                               "search",
                                               "--output search.out",
                                               "--error ",
                                               std::to_string(number_of_errors),
                                               "--index ",
                                               ibf_path(number_of_repeated_bins, window_size),
                                               "--quiet",
                                               "--query ",
                                               data("query_empty.fq"));
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);

    compare_search(number_of_repeated_bins, number_of_errors, "search.out", is_empty::yes);
}

TEST_F(search_ibf, cache_thresholds)
{
    size_t const number_of_repeated_bins{16};
    uint32_t const window_size{23};
    uint8_t const number_of_errors{1};

    {
        cli_test_result const result = execute_app("raptor",
                                                   "search",
                                                   "--cache-thresholds",
                                                   "--output search.out",
                                                   "--error ",
                                                   std::to_string(number_of_errors),
                                                   "--p_max 0.4",
                                                   "--index ",
                                                   ibf_path(number_of_repeated_bins, window_size),
                                                   "--quiet",
                                                   "--query ",
                                                   data("query.fq"));
        EXPECT_EQ(result.out, std::string{});
        EXPECT_EQ(result.err, std::string{});
        RAPTOR_ASSERT_ZERO_EXIT(result);
    }

    cli_test_result const result = execute_app("raptor",
                                               "search",
                                               "--cache-thresholds",
                                               "--output search.out",
                                               "--error ",
                                               std::to_string(number_of_errors),
                                               "--p_max 0.4",
                                               "--index ",
                                               ibf_path(number_of_repeated_bins, window_size),
                                               "--quiet",
                                               "--query ",
                                               data("query.fq"));
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);

    compare_search(number_of_repeated_bins, number_of_errors, "search.out");
}

INSTANTIATE_TEST_SUITE_P(search_ibf_suite,
                         search_ibf,
                         testing::Combine(testing::Values(0, 16, 32), testing::Values(19, 23), testing::Values(0, 1)),
                         [](testing::TestParamInfo<search_ibf::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::max<int>(1, std::get<0>(info.param) * 4)) + "_bins_"
                                              + std::to_string(std::get<1>(info.param)) + "_window_"
                                              + std::to_string(std::get<2>(info.param)) + "_error";
                             return name;
                         });

TEST_F(search_ibf, verbose)
{
    cli_test_result const result = execute_app("raptor",
                                               "search",
                                               "--output search.out",
                                               "--timing-output raptor.time",
                                               "--error 1",
                                               "--p_max 0.4",
                                               "--index ",
                                               ibf_path(16, 19),
                                               "--query ",
                                               data("query.fq"));
    EXPECT_EQ(result.out, std::string{});
    EXPECT_TRUE(result.err.starts_with("============= Timings ============="));
    EXPECT_TRUE(std::filesystem::exists("raptor.time"));
    RAPTOR_ASSERT_ZERO_EXIT(result);

    compare_search(16, 1, "search.out");
}

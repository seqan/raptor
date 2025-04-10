// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <raptor/test/cli_test.hpp>

// Same as search_ibf_test.cpp, but
// * `--fpga` as option
// * Expection result.err to be non-empty (profiling)
// * Not running w==k
// * Not running 1 bin
// * Not running threshold
// * Not running cache_tresholds
// * Not running verbose

struct search_ibf_fpga : public raptor_base, public testing::WithParamInterface<std::tuple<size_t, size_t, size_t>>
{};

TEST_P(search_ibf_fpga, with_error)
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
                                               "--fpga",
                                               "--query ",
                                               data("query.fq"));
    EXPECT_EQ(result.out, std::string{});
    EXPECT_NE(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);

    compare_search(number_of_repeated_bins, number_of_errors, "search.out");
}

TEST_P(search_ibf_fpga, no_hits)
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
                                               "--fpga",
                                               "--query ",
                                               data("query_empty.fq"));
    EXPECT_EQ(result.out, std::string{});
    EXPECT_NE(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);

    compare_search(number_of_repeated_bins, number_of_errors, "search.out", is_empty::yes);
}

INSTANTIATE_TEST_SUITE_P(search_ibf_fpga_suite,
                         search_ibf_fpga,
                         testing::Combine(testing::Values(16, 32), testing::Values(23), testing::Values(0, 1)),
                         [](testing::TestParamInfo<search_ibf_fpga::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::max<int>(1, std::get<0>(info.param) * 4)) + "_bins_"
                                              + std::to_string(std::get<1>(info.param)) + "_window_"
                                              + std::to_string(std::get<2>(info.param)) + "_error";
                             return name;
                         });

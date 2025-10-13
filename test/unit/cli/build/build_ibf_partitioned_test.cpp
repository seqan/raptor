// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h> // for Message, AssertionResult, TestParamInfo, TestPartRe...

#include <algorithm>   // for max
#include <cstddef>     // for size_t
#include <filesystem>  // for path, exists
#include <fstream>     // for operator<<, basic_ofstream, basic_ostream, ofstream
#include <ranges>      // for join_view
#include <sstream>     // for basic_stringstream
#include <string>      // for basic_string, char_traits, operator+, to_string
#include <string_view> // for basic_string_view
#include <tuple>       // for get, tuple

#include <seqan3/core/debug_stream/detail/to_string.hpp> // for to_string
#include <seqan3/test/pretty_printing.hpp>               // for PrintTo

#include <raptor/test/cli_test.hpp> // for RAPTOR_ASSERT_ZERO_EXIT, raptor_base, raptor_base::...

struct build_ibf_partitioned :
    public raptor_base,
    public testing::WithParamInterface<std::tuple<size_t, size_t, size_t, size_t>>
{};

TEST_P(build_ibf_partitioned, pipeline)
{
    auto const [number_of_repeated_bins, window_size, number_of_errors, parts] = GetParam();

    std::stringstream header{};
    { // generate input file
        std::ofstream file{"raptor_cli_test.txt"};
        size_t usr_bin_id{0};
        for (auto && file_path : get_repeated_bins(number_of_repeated_bins))
        {
            header << '#' << usr_bin_id++ << '\t' << file_path << '\n';
            file << file_path << '\n';
        }
        header << "#QUERY_NAME\tUSER_BINS\n";
        file << '\n';
    }

    cli_test_result const result1 = execute_app("raptor",
                                                "build",
                                                "--kmer 19",
                                                "--window ",
                                                std::to_string(window_size),
                                                "--output raptor.index",
                                                "--threads 2",
                                                "--parts ",
                                                std::to_string(parts),
                                                "--timing-output raptor.time",
                                                "--input",
                                                "raptor_cli_test.txt");
    EXPECT_EQ(result1.out, std::string{});
    EXPECT_TRUE(result1.err.starts_with("============= Timings ============="));
    EXPECT_TRUE(std::filesystem::exists("raptor.time"));
    RAPTOR_ASSERT_ZERO_EXIT(result1);

    cli_test_result const result2 = execute_app("raptor",
                                                "search",
                                                "--output search.out",
                                                "--error ",
                                                std::to_string(number_of_errors),
                                                "--index ",
                                                "raptor.index",
                                                "--timing-output search.time",
                                                "--query ",
                                                data("query.fq"));
    EXPECT_EQ(result2.out, std::string{});
    EXPECT_TRUE(result2.err.starts_with("============= Timings ============="));
    EXPECT_TRUE(std::filesystem::exists("search.time"));
    RAPTOR_ASSERT_ZERO_EXIT(result2);

    compare_search(number_of_repeated_bins, number_of_errors, "search.out");
}

TEST_F(build_ibf_partitioned, pipeline_misc)
{
    std::stringstream header{};
    { // generate input files
        std::ofstream file{"raptor_cli_test.txt"};
        std::ofstream file2{"raptor_cli_test.minimiser"};
        size_t usr_bin_id{0};
        for (auto && file_path : get_repeated_bins(16))
        {
            file << file_path << '\n';
            auto line = seqan3::detail::to_string("precomputed_minimisers/",
                                                  std::filesystem::path{file_path}.stem().c_str(),
                                                  ".minimiser");
            header << '#' << usr_bin_id++ << '\t' << line << '\n';
            file2 << line << '\n';
        }
        header << "#QUERY_NAME\tUSER_BINS\n";
        file << '\n';
    }

    cli_test_result const result1 = execute_app("raptor",
                                                "prepare",
                                                "--kmer 19",
                                                "--window 23",
                                                "--output precomputed_minimisers",
                                                "--quiet",
                                                "--input",
                                                "raptor_cli_test.txt");
    EXPECT_EQ(result1.out, std::string{});
    EXPECT_EQ(result1.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result1);

    cli_test_result const result2 = execute_app("raptor",
                                                "build",
                                                "--fpr 0.016",
                                                "--output raptor.index",
                                                "--parts 4",
                                                "--quiet",
                                                "--input",
                                                "raptor_cli_test.minimiser");
    EXPECT_EQ(result2.out, std::string{});
    EXPECT_EQ(result2.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result2);

    cli_test_result const result3 = execute_app("raptor",
                                                "search",
                                                "--output search.out",
                                                "--threshold 0.5",
                                                "--index ",
                                                "raptor.index",
                                                "--quiet",
                                                "--query ",
                                                data("query.fq"));
    EXPECT_EQ(result3.out, std::string{});
    EXPECT_EQ(result3.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result3);

    compare_search(16, 1 /* Always finds everything */, "search.out");

    cli_test_result const result4 = execute_app("raptor",
                                                "search",
                                                "--output search2.out",
                                                "--error 1",
                                                "--index ",
                                                "raptor.index",
                                                "--quiet",
                                                "--query ",
                                                data("query_empty.fq"));
    EXPECT_EQ(result4.out, std::string{});
    EXPECT_EQ(result4.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result4);

    compare_search(16, 1, "search2.out", is_empty::yes);
}

INSTANTIATE_TEST_SUITE_P(
    build_ibf_partitioned_suite,
    build_ibf_partitioned,
    testing::Combine(testing::Values(32), testing::Values(19, 23), testing::Values(0, 1), testing::Values(2, 4, 8)),
    [](testing::TestParamInfo<build_ibf_partitioned::ParamType> const & info)
    {
        std::string name = std::to_string(std::max<int>(1, std::get<0>(info.param) * 4)) + "_bins_"
                         + std::to_string(std::get<1>(info.param)) + "_window_"
                         + std::to_string(std::get<2>(info.param)) + "_error" + std::to_string(std::get<3>(info.param))
                         + "_parts";
        return name;
    });

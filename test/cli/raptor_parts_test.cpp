// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include "cli_test.hpp"

struct raptor_parts : public raptor_base, public testing::WithParamInterface<std::tuple<size_t, size_t, size_t, size_t, bool>> {};

TEST_P(raptor_parts, pipeline)
{
    auto const [number_of_repeated_bins, window_size, number_of_errors, parts, compressed] = GetParam();

    if (window_size == 23 && number_of_errors == 0)
        GTEST_SKIP() << "Needs dynamic threshold correction";

    std::stringstream header{};
    {
        std::string const expanded_bins = repeat_bins(number_of_repeated_bins);
        std::ofstream file{"raptor_cli_test.txt"};
        auto split_bins = expanded_bins
                        | std::views::split(' ')
                        | std::views::transform([](auto &&rng) {
                            return std::string_view(&*rng.begin(), std::ranges::distance(rng));});
        size_t usr_bin_id{0};
        for (auto && file_path : split_bins)
        {
            header << '#' << usr_bin_id++ << '\t' << file_path << '\n';
            file << file_path << '\n';
        }
        header << "#QUERY_NAME\tUSER_BINS\n";
        file << '\n';
    }

    cli_test_result const result1 = execute_app("raptor", "build",
                                                          "--kmer 19",
                                                          "--window ", std::to_string(window_size),
                                                          "--size 64k",
                                                          "--output raptor.index",
                                                          compressed ? "--compressed" : "--threads 1",
                                                          "--parts ", std::to_string(parts),
                                                          "raptor_cli_test.txt");
    EXPECT_EQ(result1.out, std::string{});
    EXPECT_EQ(result1.err, std::string{});
    ASSERT_EQ(result1.exit_code, 0);

    cli_test_result const result2 = execute_app("raptor", "search",
                                                          "--output search.out",
                                                          "--error ", std::to_string(number_of_errors),
                                                          "--index ", "raptor.index",
                                                          "--query ", data("query.fq"));
    EXPECT_EQ(result2.out, std::string{});
    EXPECT_EQ(result2.err, std::string{});
    ASSERT_EQ(result2.exit_code, 0);

    std::string const expected = [&] ()
    {
        std::string result{header.str()};
        std::string line{};
        std::ifstream search_result{search_result_path(number_of_repeated_bins,
                                                       window_size,
                                                       number_of_errors)};
        while (std::getline(search_result, line) && line.substr(0, 6) != "query1")
        {}
        result += line;
        result += '\n';
        while (std::getline(search_result, line))
        {
            result += line;
            result += '\n';
        }

        return result;
    }();

    std::string const actual = string_from_file("search.out");

    EXPECT_EQ(expected, actual);
}

TEST_F(raptor_parts, pipeline_misc)
{
    std::stringstream header{};
    {
        std::string const expanded_bins = repeat_bins(16);
        std::ofstream file{"raptor_cli_test.txt"};
        auto split_bins = expanded_bins
                        | std::views::split(' ')
                        | std::views::transform([](auto &&rng) {
                            return std::string_view(&*rng.begin(), std::ranges::distance(rng));});
        size_t usr_bin_id{0};
        for (auto && file_path : split_bins)
        {
            header << '#' << usr_bin_id++ << '\t' << file_path << '\n';
            file << file_path << '\n';
        }
        header << "#QUERY_NAME\tUSER_BINS\n";
        file << '\n';
    }

    cli_test_result const result1 = execute_app("raptor", "build",
                                                          "--kmer 19",
                                                          "--window 23",
                                                          "--size 64k",
                                                          "--output raptor.index",
                                                          "--parts 4",
                                                          "raptor_cli_test.txt");
    EXPECT_EQ(result1.out, std::string{});
    EXPECT_EQ(result1.err, std::string{});
    ASSERT_EQ(result1.exit_code, 0);

    cli_test_result const result2 = execute_app("raptor", "search",
                                                          "--output search.out",
                                                          "--threshold 0.5",
                                                          "--index ", "raptor.index",
                                                          "--query ", data("query.fq"));
    EXPECT_EQ(result2.out, std::string{});
    EXPECT_EQ(result2.err, std::string{});
    ASSERT_EQ(result2.exit_code, 0);

    std::string const expected = [&] ()
    {
        std::string const bin_list = [&] ()
        {
            std::string result;
            for (size_t i = 0; i < std::max<size_t>(1, 16 * 4u); ++i)
            {
                result += std::to_string(i);
                result += ',';
            }
            result.pop_back();
            return result;
        }();

        return header.str() + "query1\t" + bin_list + "\nquery2\t" + bin_list + "\nquery3\t" + bin_list + '\n';
    }();

    std::string const actual = string_from_file("search.out");

    EXPECT_EQ(expected, actual);

    cli_test_result const result3 = execute_app("raptor", "search",
                                                          "--output search2.out",
                                                          "--error 1",
                                                          "--index ", "raptor.index",
                                                          "--query ", data("query_empty.fq"));
    EXPECT_EQ(result3.out, std::string{});
    EXPECT_EQ(result3.err, std::string{});
    ASSERT_EQ(result3.exit_code, 0);

    std::string const expected2 = [&] ()
    {
        std::string result{header.str()};
        std::string line{};
        std::ifstream search_result{search_result_path(16,
                                                       23,
                                                       1,
                                                       false,
                                                       true)};
        while (std::getline(search_result, line) && line.substr(0, 6) != "query1")
        {}
        result += line;
        result += '\n';
        while (std::getline(search_result, line))
        {
            result += line;
            result += '\n';
        }

        return result;
    }();

    std::string const actual2 = string_from_file("search2.out");

    EXPECT_EQ(expected2, actual2);
}

INSTANTIATE_TEST_SUITE_P(parts_suite,
                         raptor_parts,
                         testing::Combine(testing::Values(32), testing::Values(19, 23), testing::Values(0, 1), testing::Values(2, 4, 8), testing::Values(true, false)),
                         [] (testing::TestParamInfo<raptor_parts::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::max<int>(1, std::get<0>(info.param) * 4)) + "_bins_" +
                                                std::to_string(std::get<1>(info.param)) + "_window_" +
                                                std::to_string(std::get<2>(info.param)) + "_error" +
                                                std::to_string(std::get<3>(info.param)) + "_parts" +
                                                (std::get<4>(info.param) ? "compressed" : "uncompressed");
                             return name;
                         });

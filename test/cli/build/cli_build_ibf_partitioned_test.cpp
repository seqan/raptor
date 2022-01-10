// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include "../cli_test.hpp"

struct build_ibf_partitioned : public raptor_base,
                               public testing::WithParamInterface<std::tuple<size_t, size_t, size_t, size_t, bool>> {};

TEST_P(build_ibf_partitioned, pipeline)
{
    auto const [number_of_repeated_bins, window_size, number_of_errors, parts, compressed] = GetParam();

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
    RAPTOR_ASSERT_ZERO_EXIT(result1);

    cli_test_result const result2 = execute_app("raptor", "search",
                                                          "--fpr 0.05",
                                                          "--output search.out",
                                                          "--error ", std::to_string(number_of_errors),
                                                          "--p_max 0.4",
                                                          "--index ", "raptor.index",
                                                          "--query ", data("query.fq"));
    EXPECT_EQ(result2.out, std::string{});
    EXPECT_EQ(result2.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result2);

    compare_search(number_of_repeated_bins, number_of_errors, "search.out");
}

TEST_F(build_ibf_partitioned, pipeline_misc)
{
    std::stringstream header{};
    { // generate input file
        std::ofstream file{"raptor_cli_test.txt"};
        size_t usr_bin_id{0};
        for (auto && file_path : get_repeated_bins(16))
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
    RAPTOR_ASSERT_ZERO_EXIT(result1);

    cli_test_result const result2 = execute_app("raptor", "search",
                                                          "--fpr 0.05",
                                                          "--output search.out",
                                                          "--threshold 0.5",
                                                          "--index ", "raptor.index",
                                                          "--query ", data("query.fq"));
    EXPECT_EQ(result2.out, std::string{});
    EXPECT_EQ(result2.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result2);

    compare_search(16, 1 /* Always finds everything */, "search.out");

    cli_test_result const result3 = execute_app("raptor", "search",
                                                          "--fpr 0.05",
                                                          "--output search2.out",
                                                          "--error 1",
                                                          "--index ", "raptor.index",
                                                          "--query ", data("query_empty.fq"));
    EXPECT_EQ(result3.out, std::string{});
    EXPECT_EQ(result3.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result3);

    compare_search(16, 1, "search2.out", is_empty::yes);
}

INSTANTIATE_TEST_SUITE_P(
    build_ibf_partitioned_suite,
    build_ibf_partitioned,
    testing::Combine(testing::Values(32),
                     testing::Values(19, 23),
                     testing::Values(0, 1),
                     testing::Values(2, 4, 8),
                     testing::Values(true, false)),
    [] (testing::TestParamInfo<build_ibf_partitioned::ParamType> const & info)
    {
        std::string name = std::to_string(std::max<int>(1, std::get<0>(info.param) * 4)) + "_bins_" +
                        std::to_string(std::get<1>(info.param)) + "_window_" +
                        std::to_string(std::get<2>(info.param)) + "_error" +
                        std::to_string(std::get<3>(info.param)) + "_parts" +
                        (std::get<4>(info.param) ? "compressed" : "uncompressed");
        return name;
    });

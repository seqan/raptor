// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <raptor/test/cli_test.hpp>

struct search_partitioned_hibf : public raptor_base, public testing::WithParamInterface<std::tuple<size_t, size_t, size_t>>
{};

TEST_P(search_partitioned_hibf, with_threshold)
{
    auto const [number_of_repeated_bins, parts, number_of_errors] = GetParam();

    std::filesystem::path const input_filename = "raptor_cli_test.txt";
    std::filesystem::path const layout_filename = "raptor_cli_test.layout";
    std::filesystem::path const index_filename = "raptor.index";
    std::filesystem::path const search_filename = "search.out";

    { // generate sequence (data) input file
        std::ofstream file{input_filename};
        size_t dummy_group{}; // to avoid clustering by filenames in chopper sketch
        for (auto && file_path : get_repeated_bins(number_of_repeated_bins))
            file << file_path << '\t' << ++dummy_group << '\n';
    }

    ASSERT_TRUE(std::filesystem::exists(input_filename));

    { // build layout
        cli_test_result const result = execute_app("raptor",
                                                   "layout",
                                                   "--kmer 19",
                                                   "--threads 1",
                                                   "--input",
                                                   input_filename,
                                                   "--tmax 64",
                                                   "--fpr 0.05",
                                                   "--relaxed-fpr 0.5",
                                                   "--number-of-partitions",
                                                   std::to_string(parts),
                                                   "--output",
                                                   layout_filename);

        EXPECT_EQ(result.out, std::string{});
        EXPECT_EQ(result.err, std::string{});
        RAPTOR_ASSERT_ZERO_EXIT(result);
    }

    ASSERT_TRUE(std::filesystem::exists(layout_filename));

    { // build index
        cli_test_result const result = execute_app("raptor",
                                                   "build",
                                                   "--output",
                                                   index_filename,
                                                   "--quiet",
                                                   "--input",
                                                   layout_filename);

        EXPECT_EQ(result.out, std::string{});
        EXPECT_EQ(result.err, std::string{});
        RAPTOR_ASSERT_ZERO_EXIT(result);
    }

    { // check with search if index contains expected input
        cli_test_result const result = execute_app("raptor",
                                                   "search",
                                                   "--output",
                                                   search_filename,
                                                   "--error ",
                                                   std::to_string(number_of_errors),
                                                   "--index",
                                                   index_filename,
                                                   "--quiet",
                                                   "--query",
                                                   data("query.fq"));
        EXPECT_EQ(result.out, std::string{});
        EXPECT_EQ(result.err, std::string{});
        RAPTOR_ASSERT_ZERO_EXIT(result);
    }

    ASSERT_TRUE(std::filesystem::exists(search_filename));
    compare_search(number_of_repeated_bins, number_of_errors, search_filename.c_str());
}

TEST_P(search_partitioned_hibf, no_hits)
{
    auto const [number_of_repeated_bins, parts, number_of_errors] = GetParam();

    std::filesystem::path const input_filename = "raptor_cli_test.txt";
    std::filesystem::path const layout_filename = "raptor_cli_test.layout";
    std::filesystem::path const index_filename = "raptor.index";
    std::filesystem::path const search_filename = "search.out";

    { // generate sequence (data) input file
        std::ofstream file{input_filename};
        size_t dummy_group{}; // to avoid clustering by filenames in chopper sketch
        for (auto && file_path : get_repeated_bins(number_of_repeated_bins))
            file << file_path << '\t' << ++dummy_group << '\n';
    }

    ASSERT_TRUE(std::filesystem::exists(input_filename));

    { // build layout
        cli_test_result const result = execute_app("raptor",
                                                   "layout",
                                                   "--kmer 19",
                                                   "--threads 1",
                                                   "--input",
                                                   input_filename,
                                                   "--tmax 64",
                                                   "--fpr 0.05",
                                                   "--relaxed-fpr 0.5",
                                                   "--number-of-partitions",
                                                   std::to_string(parts),
                                                   "--output",
                                                   layout_filename);

        EXPECT_EQ(result.out, std::string{});
        EXPECT_EQ(result.err, std::string{});
        RAPTOR_ASSERT_ZERO_EXIT(result);
    }

    ASSERT_TRUE(std::filesystem::exists(layout_filename));

    { // build index
        cli_test_result const result = execute_app("raptor",
                                                   "build",
                                                   "--output",
                                                   index_filename,
                                                   "--quiet",
                                                   "--input",
                                                   layout_filename);

        EXPECT_EQ(result.out, std::string{});
        EXPECT_EQ(result.err, std::string{});
        RAPTOR_ASSERT_ZERO_EXIT(result);
    }

    {
        cli_test_result const result = execute_app("raptor",
                                                   "search",
                                                   "--output search.out",
                                                   "--error ",
                                                   std::to_string(number_of_errors),
                                                   "--tau 0.99",
                                                   "--index ",
                                                   index_filename,
                                                   "--quiet",
                                                   "--query ",
                                                   data("query_empty.fq"));
        EXPECT_EQ(result.out, std::string{});
        EXPECT_EQ(result.err, std::string{});
        RAPTOR_ASSERT_ZERO_EXIT(result);
    }

    compare_search(number_of_repeated_bins, number_of_errors, "search.out", is_empty::yes);
}

INSTANTIATE_TEST_SUITE_P(search_partitioned_hibf_suite,
                         search_partitioned_hibf,
                         testing::Combine(testing::Values(16, 32), testing::Values(2, 4, 8), testing::Values(0, 1)),
                         [](testing::TestParamInfo<search_partitioned_hibf::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::max<int>(1, std::get<0>(info.param) * 4)) + "_bins_"
                                              + std::to_string(std::get<1>(info.param)) + "_parts"
                                              + std::to_string(std::get<2>(info.param)) + "_error";
                             return name;
                         });

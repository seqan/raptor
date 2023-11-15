// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <raptor/test/cli_test.hpp>

struct build_hibf_layout : public raptor_base
{};

TEST_F(build_hibf_layout, pipeline)
{
    std::filesystem::path const data_filename = "raptor_cli_test.txt";
    std::filesystem::path const layout_filename = "raptor_cli_test.layout";
    std::filesystem::path const index_filename = "raptor.index";
    std::filesystem::path const search_filename = "search.out";
    size_t const number_of_repeated_bins{64};
    size_t const number_of_errors{0}; // search

    { // generate sequence (data) input file
        std::ofstream file{data_filename};
        size_t dummy_group{}; // to avoid clustering by filenames in chopper sketch
        for (auto && file_path : get_repeated_bins(number_of_repeated_bins))
            file << file_path << '\t' << ++dummy_group << '\n';
    }

    ASSERT_TRUE(std::filesystem::exists(data_filename));

    { // build layout
        cli_test_result const result = execute_app("raptor",
                                                   "layout",
                                                   "--kmer 19",
                                                   "--threads 1",
                                                   "--partitioning-approach 1",
                                                   "--input",
                                                   data_filename,
                                                   "--tmax 64",
                                                   "--fpr 0.05",
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
                                                   "--threads 1",
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
                                                   "--error",
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

    compare_search(number_of_repeated_bins, number_of_errors, search_filename.c_str());
}

// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <seqan3/test/tmp_directory.hpp>

#include <chopper/count/execute.hpp>
#include <chopper/layout/execute.hpp>

#include <raptor/test/cli_test.hpp>

struct build_hibf_layout : public raptor_base
{};

TEST_F(build_hibf_layout, pipeline)
{
    seqan3::test::tmp_directory const out_dir{};
    std::filesystem::path const data_filename = out_dir.path() / "raptor_cli_test.txt";
    std::filesystem::path const layout_filename = out_dir.path() / "raptor_cli_test.layout";
    std::filesystem::path const index_filename = out_dir.path() / "raptor.index";
    std::filesystem::path const search_filename = out_dir.path() / "search.out";
    size_t const number_of_repeated_bins{16};
    size_t const number_of_errors{0}; // search

    { // generate sequence (data) input file
        std::ofstream file{data_filename};
        size_t dummy_group{}; // to avoid clustering by filenames in chopper count
        for (auto && file_path : get_repeated_bins(number_of_repeated_bins))
            file << file_path << '\t' << ++dummy_group << '\n';
    }

    ASSERT_TRUE(std::filesystem::exists(data_filename));

    { // build layout
        cli_test_result const result = execute_app("raptor",
                                                   "layout",
                                                   "--kmer-size 19",
                                                   "--column-index 2",
                                                   "--threads 1",
                                                   "--input-file",
                                                   data_filename.c_str(),
                                                   "--tmax 64",
                                                   "--false-positive-rate 0.05",
                                                   "--output-filename",
                                                   layout_filename.c_str());

        EXPECT_EQ(result.out, std::string{});
        EXPECT_EQ(result.err, std::string{});
        RAPTOR_ASSERT_ZERO_EXIT(result);
    }

    ASSERT_TRUE(std::filesystem::exists(layout_filename));

    { // build index
        cli_test_result const result = execute_app("raptor",
                                                   "build",
                                                   "--hibf",
                                                   "--kmer 19",
                                                   "--window 19",
                                                   "--fpr 0.05",
                                                   "--threads 1",
                                                   "--output",
                                                   index_filename,
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
                                                   "--query",
                                                   data("query.fq"));
        EXPECT_EQ(result.out, std::string{});
        EXPECT_EQ(result.err, std::string{});
        RAPTOR_ASSERT_ZERO_EXIT(result);
    }

    compare_search(number_of_repeated_bins, number_of_errors, search_filename.c_str());
}

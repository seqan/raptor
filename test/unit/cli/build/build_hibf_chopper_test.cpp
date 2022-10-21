// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <seqan3/test/tmp_directory.hpp>
#include <seqan3/test/tmp_filename.hpp>

#include <chopper/count/execute.hpp>
#include <chopper/layout/execute.hpp>

#include "../../../include/cli_test.hpp"

struct build_hibf_chopper : public raptor_base
{};

TEST_F(build_hibf_chopper, pipeline)
{
    seqan3::test::tmp_directory count_dir{};
    std::filesystem::path const count_prefix = count_dir.path() / "raptor_cli_test";
    seqan3::test::tmp_filename const data_filename{"raptor_cli_test.txt"};
    seqan3::test::tmp_filename const layout_filename{"raptor_cli_test.layout"};
    seqan3::test::tmp_filename const index_filename{"raptor.index"};
    seqan3::test::tmp_filename const search_filename{"search.out"};
    size_t const number_of_repeated_bins{16};
    size_t const window_size{19};
    size_t const number_of_errors{0}; // search

    { // generate sequence (data) input file
        std::ofstream file{data_filename.get_path()};
        size_t dummy_group{}; // to avoid clustering by filenames in chopper count
        for (auto && file_path : get_repeated_bins(number_of_repeated_bins))
            file << file_path << '\t' << ++dummy_group << '\n';
    }

    ASSERT_TRUE(std::filesystem::exists(data_filename.get_path()));

    { // generate count file
        char const * argv[] = {"./chopper-count",
                               "--kmer-size",
                               "19",
                               "--disable-sketch-output",
                               "--column-index",
                               "2",
                               "--threads",
                               "1",
                               "--input-file",
                               data_filename.get_path().c_str(),
                               "--output-prefix",
                               count_prefix.c_str()};
        int const argc = sizeof(argv) / sizeof(*argv);
        seqan3::argument_parser parser{"chopper-count", argc, argv, seqan3::update_notifications::off};
        chopper::count::execute(parser);
    }

    ASSERT_TRUE(std::filesystem::exists(count_prefix.string() + ".count"));

    { // generate layout file
        char const * argv[] = {"./chopper-layout",
                               "--tmax",
                               "64",
                               "--false-positive-rate",
                               "0.05",
                               "--input-prefix",
                               count_prefix.c_str(),
                               "--output-filename",
                               layout_filename.get_path().c_str()};
        int const argc = sizeof(argv) / sizeof(*argv);
        seqan3::argument_parser parser{"chopper-layout", argc, argv, seqan3::update_notifications::off};
        chopper::layout::execute(parser);
    }

    ASSERT_TRUE(std::filesystem::exists(layout_filename.get_path()));

    { // build index
        cli_test_result const result = execute_app("raptor",
                                                   "build",
                                                   "--hibf",
                                                   "--kmer 19",
                                                   "--window",
                                                   std::to_string(window_size),
                                                   "--fpr 0.05",
                                                   "--threads 1",
                                                   "--output",
                                                   index_filename.get_path(),
                                                   layout_filename.get_path());

        EXPECT_EQ(result.out, std::string{});
        EXPECT_EQ(result.err, std::string{});
        RAPTOR_ASSERT_ZERO_EXIT(result);
    }

    { // check with search if index contains expected input
        cli_test_result const result = execute_app("raptor",
                                                   "search",
                                                   "--fpr 0.05",
                                                   "--output",
                                                   search_filename.get_path(),
                                                   "--error",
                                                   std::to_string(number_of_errors),
                                                   "--hibf",
                                                   "--index",
                                                   index_filename.get_path(),
                                                   "--query",
                                                   data("query.fq"));
        EXPECT_EQ(result.out, std::string{});
        EXPECT_EQ(result.err, std::string{});
        RAPTOR_ASSERT_ZERO_EXIT(result);
    }

    compare_search(number_of_repeated_bins, number_of_errors, search_filename.get_path().c_str());
    count_dir.clean();
}

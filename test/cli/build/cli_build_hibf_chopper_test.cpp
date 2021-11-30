// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <chopper/count/execute.hpp>
#include <chopper/layout/execute.hpp>

#include <seqan3/test/tmp_filename.hpp>

#include "../cli_test.hpp"

struct build_hibf_chopper : public raptor_base {};

TEST_F(build_hibf_chopper, pipeline)
{
    seqan3::test::tmp_filename const data_filename{"raptor_cli_test.txt"};
    seqan3::test::tmp_filename const count_filename{"raptor_cli_test.counts"};
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
        const char * argv[] = {"./chopper-count",
                               "--kmer-size", "19",
                               "--window-size", "19",
                               "--column-index", "2",
                               "--threads", "1",
                               "--data_file", data_filename.get_path().c_str(),
                               "-o", count_filename.get_path().c_str()};
        int const argc = sizeof(argv) / sizeof(*argv);
        seqan3::argument_parser parser{"chopper-count", argc, argv, seqan3::update_notifications::off};
        chopper::count::execute(parser);
    }

    ASSERT_TRUE(std::filesystem::exists(count_filename.get_path()));

    { // generate layout file
        const char * argv[] = {"./chopper-layout",
                               "--technical-bins", "64",
                               "--false-positive-rate", "0.05",
                               "--filenames", count_filename.get_path().c_str(),
                               "-o", layout_filename.get_path().c_str()};
        int const argc = sizeof(argv) / sizeof(*argv);
        seqan3::argument_parser parser{"chopper-layout", argc, argv, seqan3::update_notifications::off};
        chopper::layout::execute(parser);
    }

    ASSERT_TRUE(std::filesystem::exists(layout_filename.get_path()));

    { // build index
        cli_test_result const result = execute_app("raptor", "build",
                                                             "--hibf",
                                                             "--kmer 19",
                                                             "--window", std::to_string(window_size),
                                                             "--fpr 0.05",
                                                             "--threads 1",
                                                             "--output", index_filename.get_path(),
                                                             layout_filename.get_path());

        EXPECT_EQ(result.out, std::string{});
        EXPECT_EQ(result.err, std::string{});
        RAPTOR_ASSERT_ZERO_EXIT(result);
    }

    { // check with search if index contains expected input
        cli_test_result const result = execute_app("raptor", "search",
                                                             "--output", search_filename.get_path(),
                                                             "--error", std::to_string(number_of_errors),
                                                             "--hibf",
                                                             "--index", index_filename.get_path(),
                                                             "--query", data("query.fq"));
        EXPECT_EQ(result.out, std::string{});
        EXPECT_EQ(result.err, std::string{});
        RAPTOR_ASSERT_ZERO_EXIT(result);
    }

    compare_search(number_of_repeated_bins, number_of_errors, search_filename.get_path().c_str());
}

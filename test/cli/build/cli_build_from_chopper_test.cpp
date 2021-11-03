// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <cstring>

#include "../cli_test.hpp"

#include <chopper/count/execute.hpp>
#include <chopper/layout/execute.hpp>

#include <seqan3/test/tmp_filename.hpp>

struct build_from_chopper : public raptor_base {};

TEST_F(build_from_chopper, pipeline)
{
    seqan3::test::tmp_filename data_filename{"raptor_cli_test.txt"};
    seqan3::test::tmp_filename count_filename{"raptor_cli_test.counts"};
    seqan3::test::tmp_filename layout_filename{"raptor_cli_test.layout"};

    {// generate sequence (data) input file
        std::string const expanded_bins = repeat_bins(16);
        std::ofstream file{data_filename.get_path()};
        auto split_bins = expanded_bins
                        | std::views::split(' ')
                        | std::views::transform([](auto &&rng) {
                            return std::string_view(&*rng.begin(), std::ranges::distance(rng));});
        size_t i{0}; // dummy spec
        for (auto && file_path : split_bins)
            file << file_path << '\t' << i <<'\n';
        file << '\n';
    }

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

    cli_test_result const result = execute_app("raptor", "build",
                                                         "--hibf",
                                                         "--kmer 19",
                                                         "--window 19",
                                                         "--fpr 0.05",
                                                         "--threads 1",
                                                         "--output raptor.index",
                                                         layout_filename.get_path().string());

    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
    ASSERT_EQ(result.exit_code, 0);

    compare_results<raptor::index_structure::hibf>(data("64bins19window.index"), "raptor.index");
}

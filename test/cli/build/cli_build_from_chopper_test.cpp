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
    seqan3::test::tmp_filename index_filename{"raptor.index"};
    seqan3::test::tmp_filename search_filename{"search.out"};
    size_t number_of_repeated_bins{16};
    size_t window_size{19};
    size_t number_of_errors{0}; // search

    {// generate sequence (data) input file
        std::ofstream file{data_filename.get_path()};
        size_t i{0}; // dummy spec
        for (auto && file_path : get_repeated_bins(number_of_repeated_bins))
            file << file_path << '\t' << ++i <<'\n';
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
                                                            "--output ", index_filename.get_path().c_str(),
                                                            layout_filename.get_path().c_str());

        EXPECT_EQ(result.out, std::string{});
        EXPECT_EQ(result.err, std::string{});
        ASSERT_EQ(result.exit_code, 0);
    }

    { // check with search if index contains expected input
        cli_test_result const result = execute_app("raptor", "search",
                                                            "--output", search_filename.get_path().c_str(),
                                                            "--error ", std::to_string(number_of_errors),
                                                            "--hibf",
                                                            "--index ", index_filename.get_path().c_str(),
                                                            "--query ", data("query.fq"));
        EXPECT_EQ(result.exit_code, 0);
        EXPECT_EQ(result.out, std::string{});
        EXPECT_EQ(result.err, std::string{});

        std::ifstream in{search_filename.get_path()};
        std::string line;
        std::string query_result;

        while (std::getline(in, line) && line[0] == '#' && line[1] != 'Q')
            if (line.substr(line.size() - 4)[0] != '4') // only file bin4.fa does not contain the query
                query_result += line.substr(1, std::find(line.begin(), line.end(), '\t') - line.begin() - 1) + ",";
        query_result.pop_back(); // remove last '.'

        EXPECT_TRUE(std::getline(in, line)); // get query1
        EXPECT_EQ(line.substr(7), query_result);
        EXPECT_TRUE(std::getline(in, line)); // get query2
        EXPECT_EQ(line.substr(7), query_result);
        EXPECT_TRUE(std::getline(in, line)); // get query3
        EXPECT_EQ(line.substr(7), query_result);

        EXPECT_FALSE(std::getline(in, line));
    }
}

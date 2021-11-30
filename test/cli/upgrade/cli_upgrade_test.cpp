// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include "../cli_test.hpp"

struct upgrade : public raptor_base {};

TEST_F(upgrade, ibf)
{
    { // generate input file
        std::ofstream file{"raptor_cli_test.txt"};
        for (auto && file_path : get_repeated_bins(16))
        {
            file << file_path << '\n';
        }
        file << '\n';
    }

    cli_test_result const result = execute_app("raptor", "upgrade",
                                                         "--kmer 19",
                                                         "--window 23",
                                                         "--bins raptor_cli_test.txt",
                                                         "--input ", data("1_1.index"),
                                                         "--output raptor.index");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_RESULT(result);

    compare_results(ibf_path(16, 23), "raptor.index");
}

TEST_F(upgrade, compressed_ibf)
{
    { // generate input file
        std::ofstream file{"raptor_cli_test.txt"};
        for (auto && file_path : get_repeated_bins(16))
        {
            file << file_path << '\n';
        }
        file << '\n';
    }

    cli_test_result const result = execute_app("raptor", "upgrade",
                                                         "--kmer 19",
                                                         "--window 23",
                                                         "--compressed",
                                                         "--bins raptor_cli_test.txt",
                                                         "--input ", data("1_1c.index"),
                                                         "--output raptor.index");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_RESULT(result);

    compare_results<raptor::index_structure::ibf_compressed>(ibf_path(16, 23, true), "raptor.index");
}

TEST_F(upgrade, partitioned_ibf)
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

    cli_test_result const result1 = execute_app("raptor", "upgrade",
                                                          "--kmer 19",
                                                          "--window 23",
                                                          "--bins raptor_cli_test.txt",
                                                          "--input ", data("1_1.index"),
                                                          "--output raptor.index",
                                                          "--parts 4");

    EXPECT_EQ(result1.out, std::string{});
    EXPECT_EQ(result1.err, std::string{});
    ASSERT_EQ(result1.exit_code, 0);

    cli_test_result const result2 = execute_app("raptor", "search",
                                                          "--output search.out",
                                                          "--error 1",
                                                          "--index raptor.index",
                                                          "--query ", data("query.fq"));
    EXPECT_EQ(result2.out, std::string{});
    EXPECT_EQ(result2.err, std::string{});
    RAPTOR_ASSERT_RESULT(result2);

    std::string const expected = [&] ()
    {
        std::string result{header.str()};
        std::string line{};
        std::ifstream search_result{search_result_path(16, 23, 1)};
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

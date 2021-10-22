// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include "../cli_test.hpp"

struct build_ibf_compressed : public raptor_base, public testing::WithParamInterface<std::tuple<size_t, size_t, bool, size_t>> {};

TEST_P(build_ibf_compressed, pipeline)
{
    auto const [number_of_repeated_bins, window_size, run_parallel_tmp, number_of_errors] = GetParam();
    bool const run_parallel = run_parallel_tmp && number_of_repeated_bins >= 32;

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
                                                          "--threads ", run_parallel ? "2" : "1",
                                                          "--output raptor.index",
                                                          "--compressed",
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

TEST_P(build_ibf_compressed, pipeline_socks)
{
    auto const [number_of_repeated_bins, window_size, run_parallel_tmp, number_of_errors] = GetParam();
    bool const run_parallel = run_parallel_tmp && number_of_repeated_bins >= 32;

    if (window_size == 23 || number_of_errors != 0)
        GTEST_SKIP() << "SOCKS only supports exact kmers";

    {
        std::string const expanded_bins = repeat_bins(number_of_repeated_bins);
        std::ofstream file{"raptor_cli_test.txt"};
        auto split_bins = expanded_bins
                        | std::views::split(' ')
                        | std::views::transform([](auto &&rng) {
                            return std::string_view(&*rng.begin(), std::ranges::distance(rng));});
        for (auto && file_path : split_bins)
        {
            file << "dummy_color: " << file_path << '\n';
        }
        file << '\n';
    }

    cli_test_result const result1 = execute_app("raptor", "socks", "build",
                                                          "--kmer 19",
                                                          "--window ", std::to_string(window_size),
                                                          "--size 64k",
                                                          "--threads ", run_parallel ? "2" : "1",
                                                          "--output raptor.index",
                                                          "--compressed",
                                                          "raptor_cli_test.txt");
    EXPECT_EQ(result1.out, std::string{});
    EXPECT_EQ(result1.err, std::string{});
    ASSERT_EQ(result1.exit_code, 0);

    cli_test_result const result2 = execute_app("raptor", "socks", "lookup-kmer",
                                                          "--output search.out",
                                                          "--index ", "raptor.index",
                                                          "--query ", data("query_socks.fq"));
    EXPECT_EQ(result2.out, std::string{});
    EXPECT_EQ(result2.err, std::string{});
    ASSERT_EQ(result2.exit_code, 0);

    std::string const expected = string_from_file(search_result_path(number_of_repeated_bins,
                                                                     window_size,
                                                                     number_of_errors,
                                                                     true),
                                                  std::ios::binary);
    std::string const actual = string_from_file("search.out");

    EXPECT_EQ(expected, actual);
}

INSTANTIATE_TEST_SUITE_P(
    build_ibf_compressed_suite,
    build_ibf_compressed,
    testing::Combine(testing::Values(0, 16, 32), testing::Values(19, 23), testing::Values(true, false), testing::Values(0, 1)),
    [] (testing::TestParamInfo<build_ibf_compressed::ParamType> const & info)
    {
        std::string name = std::to_string(std::max<int>(1, std::get<0>(info.param) * 4)) + "_bins_" +
                        std::to_string(std::get<1>(info.param)) + "_window_" +
                        (std::get<2>(info.param) ? "parallel" : "serial") +
                        std::to_string(std::get<3>(info.param)) + "_error";
        return name;
    });

TEST_F(build_ibf_compressed, wrong_compression)
{
    raptor::raptor_index<raptor::index_structure::ibf_compressed> index{};

    std::ifstream is{data("1bins19window.index"), std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    EXPECT_THROW(iarchive(index), seqan3::argument_parser_error);
}

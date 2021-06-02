#include <fstream>
#include <seqan3/std/ranges>     // range comparisons
#include <string>                // strings
#include <vector>                // vectors

#include "cli_test.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// raptor build tests ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEST_P(raptor_build, build_with_file)
{
    auto const [number_of_repeated_bins, window_size, run_parallel_tmp] = GetParam();
    bool const run_parallel = run_parallel_tmp && number_of_repeated_bins >= 32;

    {
        std::string const expanded_bins = repeat_bins(number_of_repeated_bins);
        std::ofstream file{"raptor_cli_test.txt"};
        auto split_bins = expanded_bins
                        | std::views::split(' ')
                        | std::views::transform([](auto &&rng) {
                            return std::string_view(&*rng.begin(), std::ranges::distance(rng));});
        for (auto && file_path : split_bins)
        {
            file << file_path << '\n';
        }
        file << '\n';
    }

    cli_test_result const result = execute_app("raptor", "build",
                                                         "--kmer 19",
                                                         "--window ", std::to_string(window_size),
                                                         "--size 64k",
                                                         "--threads ", run_parallel ? "2" : "1",
                                                         "--output index.ibf",
                                                         "raptor_cli_test.txt");
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});

    std::string const expected = string_from_file(ibf_path(number_of_repeated_bins, window_size), std::ios::binary);
    std::string const actual = string_from_file("index.ibf", std::ios::binary);

    EXPECT_TRUE(expected == actual);
}

TEST_P(raptor_build, build_with_file_socks)
{
    auto const [number_of_repeated_bins, window_size, run_parallel_tmp] = GetParam();
    bool const run_parallel = run_parallel_tmp && number_of_repeated_bins >= 32;

    {
        std::string const expanded_bins = repeat_bins(number_of_repeated_bins);
        std::ofstream file{"raptor_cli_test.txt"};
        auto split_bins = expanded_bins
                        | std::views::split(' ')
                        | std::views::transform([](auto &&rng) {
                            return std::string_view(&*rng.begin(), std::ranges::distance(rng));});
        bool repeat{false};
        for (auto && file_path : split_bins)
        {
            file << "dummy_color: " << file_path;
            if (repeat)
                file << ' ' << file_path;
            file << '\n';
            repeat = !repeat;
        }
        file << '\n';
    }

    cli_test_result const result = execute_app("raptor", "socks", "build",
                                                         "--kmer 19",
                                                         "--window ", std::to_string(window_size),
                                                         "--size 64k",
                                                         "--threads ", run_parallel ? "2" : "1",
                                                         "--output index.ibf",
                                                         "raptor_cli_test.txt");
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});

    std::string const expected = string_from_file(ibf_path(number_of_repeated_bins, window_size), std::ios::binary);
    std::string const actual = string_from_file("index.ibf", std::ios::binary);

    EXPECT_TRUE(expected == actual);
}

INSTANTIATE_TEST_SUITE_P(build_suite,
                         raptor_build,
                         testing::Combine(testing::Values(0, 16, 32), testing::Values(19, 23), testing::Values(true, false)),
                         [] (testing::TestParamInfo<raptor_build::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::max<int>(1, std::get<0>(info.param) * 4)) + "_bins_" +
                                                std::to_string(std::get<1>(info.param)) + "_window_" +
                                                (std::get<2>(info.param) ? "parallel" : "serial");
                             return name;
                         });

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// raptor search tests //////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEST_P(raptor_search, search)
{
    auto const [number_of_repeated_bins, window_size, number_of_errors] = GetParam();

    if (window_size == 23 && number_of_errors == 0)
        GTEST_SKIP() << "Needs dynamic threshold correction";

    cli_test_result const result = execute_app("raptor", "search",
                                                         "--output search.out",
                                                         "--error ", std::to_string(number_of_errors),
                                                         "--index ", ibf_path(number_of_repeated_bins, window_size),
                                                         "--query ", data("query.fq"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});

    std::string const expected = string_from_file(search_result_path(number_of_repeated_bins, window_size, number_of_errors), std::ios::binary);
    std::string const actual = string_from_file("search.out");

    EXPECT_EQ(expected, actual);
}

TEST_P(raptor_search, search_socks)
{
    auto const [number_of_repeated_bins, window_size, number_of_errors] = GetParam();

    if (window_size == 23 || number_of_errors != 0)
        GTEST_SKIP() << "SOCKS only supports exact kmers";

    cli_test_result const result = execute_app("raptor", "socks", "lookup-kmer",
                                                         "--output search.out",
                                                         "--index ", ibf_path(number_of_repeated_bins, window_size),
                                                         "--query ", data("query_socks.fq"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});

    std::string const expected = string_from_file(search_result_path(number_of_repeated_bins, window_size, number_of_errors, true), std::ios::binary);
    std::string const actual = string_from_file("search.out");

    EXPECT_EQ(expected, actual);
}

// Search with threshold
TEST_P(raptor_search, search_threshold)
{
    auto const [number_of_repeated_bins, window_size, number_of_errors] = GetParam();

    cli_test_result const result = execute_app("raptor", "search",
                                                         "--output search_threshold.out",
                                                         "--threshold 0.50",
                                                         "--index ", ibf_path(number_of_repeated_bins, window_size),
                                                         "--query ", data("query.fq"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});

    std::string const expected = [&] ()
    {
        std::string const bin_list = [&] ()
        {
            std::string result;
            for (size_t i = 0; i < std::max<size_t>(1, number_of_repeated_bins * 4u); ++i)
            {
                result += std::to_string(i);
                result += ',';
            }
            return result;
        }();

        return "query1\t" + bin_list + "\nquery2\t" + bin_list + "\nquery3\t" + bin_list + '\n';
    }();

    std::string const actual = string_from_file("search_threshold.out");

    EXPECT_EQ(expected, actual);
}

INSTANTIATE_TEST_SUITE_P(search_suite,
                         raptor_search,
                         testing::Combine(testing::Values(0, 16, 32), testing::Values(19, 23), testing::Values(0, 1)),
                         [] (testing::TestParamInfo<raptor_search::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::max<int>(1, std::get<0>(info.param) * 4)) + "_bins_" +
                                                std::to_string(std::get<1>(info.param)) + "_window_" +
                                                std::to_string(std::get<2>(info.param)) + "_error";
                             return name;
                         });

#include <fstream>
#include <seqan3/std/ranges>     // range comparisons
#include <string>                // strings
#include <vector>                // vectors

#include "cli_test.hpp"

inline std::string string_from_file(std::filesystem::path const & path,
                                    std::ios_base::openmode const mode = std::ios_base::in)
{
    std::ifstream file_stream(path, mode);
    if (!file_stream.is_open())
        throw std::logic_error{"Cannot open " + path.string()};
    std::stringstream file_buffer;
    file_buffer << file_stream.rdbuf();
    return {file_buffer.str()};
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// raptor build tests ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEST_P(raptor_build, build_with_list)
{
    auto const [number_of_bins, run_parallel] = GetParam();

    cli_test_result const result = execute_app("raptor", "build",
                                                         "--kmer 19",
                                                         "--window 23",
                                                         "--size 8m",
                                                         "--threads ", run_parallel ? "2" : "1",
                                                         "--output index.ibf", expand_bins(number_of_bins));
    ASSERT_EQ(result.exit_code, 0);
    ASSERT_EQ(result.out, std::string{});
    ASSERT_EQ(result.err, std::string{});

    std::string const expected = string_from_file(data("expected_results/b" +
                                                       std::to_string(number_of_bins) +
                                                       "_k19_w23_s8m.ibf"),
                                                  std::ios::binary);
    std::string const actual = string_from_file("index.ibf", std::ios::binary);

    EXPECT_TRUE(expected == actual);
}

TEST_P(raptor_build, build_with_file)
{
    auto const [number_of_bins, run_parallel] = GetParam();

    {
        std::string const expanded_bins = expand_bins(number_of_bins);
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
                                                         "--window 23",
                                                         "--size 8m",
                                                         "--threads ", run_parallel ? "2" : "1",
                                                         "--output index.ibf",
                                                         "raptor_cli_test.txt");
    ASSERT_EQ(result.exit_code, 0);
    ASSERT_EQ(result.out, std::string{});
    ASSERT_EQ(result.err, std::string{});

    std::string const expected = string_from_file(data("expected_results/b" +
                                                       std::to_string(number_of_bins) +
                                                       "_k19_w23_s8m.ibf"),
                                                  std::ios::binary);
    std::string const actual = string_from_file("index.ibf", std::ios::binary);

    EXPECT_TRUE(expected == actual);
}

INSTANTIATE_TEST_SUITE_P(build_suite,
                         raptor_build,
                         testing::Combine(testing::Values(64, 1024), testing::Values(true, false)),
                         [] (testing::TestParamInfo<raptor_build::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::get<0>(info.param)) + "bins_" +
                                                (std::get<1>(info.param) ? "parallel" : "serial");
                             return name;
                         });

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// raptor search tests //////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEST_P(raptor_search, search)
{
    auto const [number_of_bins, window_size, number_of_errors] = GetParam();

    cli_test_result const result = execute_app("raptor", "search",
                                                         "--kmer 19",
                                                         "--window ", std::to_string(window_size),
                                                         "--output search.out",
                                                         "--error ", std::to_string(number_of_errors),
                                                         "--pattern 100",
                                                         "--index ", data("expected_results/b" +
                                                                          std::to_string(number_of_bins) +
                                                                          "_k19_w" +
                                                                          std::to_string(window_size) +
                                                                          "_s8m.ibf"),
                                                         "--query ", data("example_data/" +
                                                                          std::to_string(number_of_bins) +
                                                                          "/reads/all.fastq"));
    ASSERT_EQ(result.exit_code, 0);
    ASSERT_EQ(result.out, std::string{});
    ASSERT_EQ(result.err, std::string{});

    std::string const expected = string_from_file(data("expected_results/b" +
                                                       std::to_string(number_of_bins) +
                                                       "_k19_w" +
                                                       std::to_string(window_size) +
                                                       "_s8m_e" +
                                                       std::to_string(number_of_errors) +
                                                       ".out"));
    std::string const actual = string_from_file("search.out");

    EXPECT_TRUE(expected == actual);
}

// Search with threshold
TEST_P(raptor_search, search_threshold)
{
    auto const [number_of_bins, window_size, number_of_errors] = GetParam();

    cli_test_result const result = execute_app("raptor", "search",
                                                        "--kmer 19",
                                                        "--window ", std::to_string(window_size),
                                                        "--output search_threshold.out",
                                                        "--error ", std::to_string(number_of_errors),
                                                        "--threshold 0.50",
                                                        "--pattern 100",
                                                        "--index ", data("expected_results/b" +
                                                                         std::to_string(number_of_bins) +
                                                                         "_k19_w" +
                                                                         std::to_string(window_size) +
                                                                         "_s8m.ibf"),
                                                        "--query ", data("example_data/" +
                                                                         std::to_string(number_of_bins) +
                                                                         "/reads/all.fastq"));
    ASSERT_EQ(result.exit_code, 0);
    ASSERT_EQ(result.out, std::string{});
    ASSERT_EQ(result.err, std::string{});

    std::string const expected = string_from_file(data("expected_results/b" +
                                                       std::to_string(number_of_bins) +
                                                       "_k19_w" +
                                                       std::to_string(window_size) +
                                                       "_s8m_e" +
                                                       std::to_string(number_of_errors) +
                                                       ".out"));
    std::string const actual = string_from_file("search_threshold.out");

    // EXPECT_TRUE(expected == actual);
}

INSTANTIATE_TEST_SUITE_P(search_suite,
                         raptor_search,
                         testing::Combine(testing::Values(64, 1024), testing::Values(19, 23), testing::Values(0, 1, 2, 3)),
                         [] (testing::TestParamInfo<raptor_search::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::get<0>(info.param)) + "bins_" +
                                                std::to_string(std::get<1>(info.param)) + "window_" +
                                                std::to_string(std::get<2>(info.param)) + "error";
                             return name;
                         });

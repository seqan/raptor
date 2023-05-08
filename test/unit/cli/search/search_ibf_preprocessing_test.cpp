// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <raptor/test/cli_test.hpp>

struct search_ibf_preprocessing :
    public raptor_base,
    public testing::WithParamInterface<std::tuple<size_t, size_t, bool, size_t>>
{};

TEST_P(search_ibf_preprocessing, pipeline)
{
    auto const [number_of_repeated_bins, window_size, run_parallel_tmp, number_of_errors] = GetParam();
    bool const run_parallel = run_parallel_tmp && number_of_repeated_bins >= 32;

    std::stringstream header{};
    { // generate input files
        std::ofstream file{"raptor_cli_test.txt"};
        std::ofstream file2{"raptor_cli_test.minimiser"};
        size_t usr_bin_id{0};
        for (auto && file_path : get_repeated_bins(number_of_repeated_bins))
        {
            file << file_path << '\n';
            auto line = seqan3::detail::to_string("precomputed_minimisers/",
                                                  std::filesystem::path{file_path}.stem().c_str(),
                                                  ".minimiser");
            header << '#' << usr_bin_id++ << '\t' << line << '\n';
            file2 << line << '\n';
        }
        header << "#QUERY_NAME\tUSER_BINS\n";
        file << '\n';
    }

    cli_test_result const result1 = execute_app("raptor",
                                                "prepare",
                                                "--kmer 19",
                                                "--window ",
                                                std::to_string(window_size),
                                                "--threads ",
                                                run_parallel ? "2" : "1",
                                                "--output precomputed_minimisers",
                                                "--quiet",
                                                "--input raptor_cli_test.txt");
    EXPECT_EQ(result1.out, std::string{});
    EXPECT_EQ(result1.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result1);

    cli_test_result const result2 = execute_app("raptor",
                                                "build",
                                                "--threads ",
                                                run_parallel ? "2" : "1",
                                                "--output raptor.index",
                                                "--quiet",
                                                "--input raptor_cli_test.minimiser");
    EXPECT_EQ(result2.out, std::string{});
    EXPECT_EQ(result2.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result2);

    compare_index(ibf_path(number_of_repeated_bins, window_size), "raptor.index", compare_extension::no);

    cli_test_result const result3 = execute_app("raptor",
                                                "search",
                                                "--output search.out",
                                                "--error ",
                                                std::to_string(number_of_errors),
                                                "--p_max 0.4",
                                                "--index ",
                                                "raptor.index",
                                                "--quiet",
                                                "--query ",
                                                data("query.fq"));
    EXPECT_EQ(result3.out, std::string{});
    EXPECT_EQ(result3.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result3);

    compare_search(number_of_repeated_bins, number_of_errors, "search.out", is_empty::no, is_preprocessed::yes);
}

TEST_P(search_ibf_preprocessing, pipeline_compressed_bins)
{
    auto const [number_of_repeated_bins, window_size, run_parallel_tmp, number_of_errors] = GetParam();
    bool const run_parallel = run_parallel_tmp && number_of_repeated_bins >= 32;

    std::stringstream header{};
    { // generate input files
        std::ofstream file{"raptor_cli_test.txt"};
        std::ofstream file2{"raptor_cli_test.minimiser"};
        size_t usr_bin_id{0};
        for (auto && file_path : get_repeated_bins(number_of_repeated_bins))
        {
            file << file_path << ".gz\n";
            auto line = seqan3::detail::to_string("precomputed_minimisers/",
                                                  std::filesystem::path{file_path}.stem().c_str(),
                                                  ".minimiser");
            header << '#' << usr_bin_id++ << '\t' << line << '\n';
            file2 << line << '\n';
        }
        header << "#QUERY_NAME\tUSER_BINS\n";
        file << '\n';
    }

    cli_test_result const result1 = execute_app("raptor",
                                                "prepare",
                                                "--kmer 19",
                                                "--window ",
                                                std::to_string(window_size),
                                                "--use-filesize-dependent-cutoff",
                                                "--threads ",
                                                run_parallel ? "2" : "1",
                                                "--output precomputed_minimisers",
                                                "--quiet",
                                                "--input raptor_cli_test.txt");
    EXPECT_EQ(result1.out, std::string{});
    EXPECT_EQ(result1.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result1);

    cli_test_result const result2 = execute_app("raptor",
                                                "build",
                                                "--threads ",
                                                run_parallel ? "2" : "1",
                                                "--output raptor.index",
                                                "--quiet",
                                                "--input raptor_cli_test.minimiser");
    EXPECT_EQ(result2.out, std::string{});
    EXPECT_EQ(result2.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result2);

    compare_index(ibf_path(number_of_repeated_bins, window_size), "raptor.index", compare_extension::no);

    cli_test_result const result3 = execute_app("raptor",
                                                "search",
                                                "--output search.out",
                                                "--error ",
                                                std::to_string(number_of_errors),
                                                "--index ",
                                                "raptor.index",
                                                "--quiet",
                                                "--query ",
                                                data("query.fq"));
    EXPECT_EQ(result3.out, std::string{});
    EXPECT_EQ(result3.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result3);

    compare_search(number_of_repeated_bins, number_of_errors, "search.out", is_empty::no, is_preprocessed::yes);
}

INSTANTIATE_TEST_SUITE_P(search_ibf_preprocessing_suite,
                         search_ibf_preprocessing,
                         testing::Combine(testing::Values(0, 16, 32),
                                          testing::Values(19, 23),
                                          testing::Values(true, false),
                                          testing::Values(0, 1)),
                         [](testing::TestParamInfo<search_ibf_preprocessing::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::max<int>(1, std::get<0>(info.param) * 4)) + "_bins_"
                                              + std::to_string(std::get<1>(info.param)) + "_window_"
                                              + (std::get<2>(info.param) ? "parallel" : "serial")
                                              + std::to_string(std::get<3>(info.param)) + "_error";
                             return name;
                         });

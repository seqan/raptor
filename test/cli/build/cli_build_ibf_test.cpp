// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include "../cli_test.hpp"

struct build_ibf : public raptor_base, public testing::WithParamInterface<std::tuple<size_t, size_t, bool>>
{};

TEST_P(build_ibf, with_file)
{
    auto const [number_of_repeated_bins, window_size, run_parallel_tmp] = GetParam();
    bool const run_parallel = run_parallel_tmp && number_of_repeated_bins >= 32;

    { // generate input file
        std::ofstream file{"raptor_cli_test.txt"};
        for (auto && file_path : get_repeated_bins(number_of_repeated_bins))
            file << file_path << '\n';
        file << '\n';
    }

    cli_test_result const result = execute_app("raptor",
                                               "build",
                                               "--kmer 19",
                                               "--window ",
                                               std::to_string(window_size),
                                               "--size 64k",
                                               "--threads ",
                                               run_parallel ? "2" : "1",
                                               "--output raptor.index",
                                               "raptor_cli_test.txt");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);

    compare_index(ibf_path(number_of_repeated_bins, window_size), "raptor.index");
}

TEST_P(build_ibf, with_shape)
{
    auto const [number_of_repeated_bins, window_size, run_parallel_tmp] = GetParam();
    bool const run_parallel = run_parallel_tmp && number_of_repeated_bins >= 32;

    { // generate input file
        std::ofstream file{"raptor_cli_test.txt"};
        for (auto && file_path : get_repeated_bins(number_of_repeated_bins))
            file << file_path << '\n';
        file << '\n';
    }

    cli_test_result const result = execute_app("raptor",
                                               "build",
                                               "--shape 1111111111111111111",
                                               "--window ",
                                               std::to_string(window_size),
                                               "--size 64k",
                                               "--threads ",
                                               run_parallel ? "2" : "1",
                                               "--output raptor.index",
                                               "raptor_cli_test.txt");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);

    compare_index(ibf_path(number_of_repeated_bins, window_size), "raptor.index");
}

TEST_P(build_ibf, with_socks_file)
{
    auto const [number_of_repeated_bins, window_size, run_parallel_tmp] = GetParam();
    bool const run_parallel = run_parallel_tmp && number_of_repeated_bins >= 32;

    { // generate input file
        std::ofstream file{"raptor_cli_test.txt"};
        bool repeat{false};
        for (auto && file_path : get_repeated_bins(number_of_repeated_bins))
        {
            file << "dummy_color: " << file_path;
            if (repeat)
                file << ' ' << file_path;
            file << '\n';
            repeat = !repeat;
        }
        file << '\n';
    }

    cli_test_result const result = execute_app("raptor",
                                               "socks",
                                               "build",
                                               "--kmer 19",
                                               "--window ",
                                               std::to_string(window_size),
                                               "--size 64k",
                                               "--threads ",
                                               run_parallel ? "2" : "1",
                                               "--output raptor.index",
                                               "raptor_cli_test.txt");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);

    compare_index(ibf_path(number_of_repeated_bins, window_size), "raptor.index");
}

INSTANTIATE_TEST_SUITE_P(build_ibf_suite,
                         build_ibf,
                         testing::Combine(testing::Values(0, 16, 32),
                                          testing::Values(19, 23),
                                          testing::Values(true, false)),
                         [](testing::TestParamInfo<build_ibf::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::max<int>(1, std::get<0>(info.param) * 4)) + "_bins_"
                                              + std::to_string(std::get<1>(info.param)) + "_window_"
                                              + (std::get<2>(info.param) ? "parallel" : "serial");
                             return name;
                         });

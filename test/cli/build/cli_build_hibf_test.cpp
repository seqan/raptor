// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include "../cli_test.hpp"

struct build_hibf : public raptor_base, public testing::WithParamInterface<std::tuple<size_t, size_t, bool>> {};

TEST_P(build_hibf, with_file)
{
    auto const [number_of_repeated_bins, window_size, run_parallel_tmp] = GetParam();
    bool const run_parallel = run_parallel_tmp && number_of_repeated_bins >= 32;

    cli_test_result const result = execute_app("raptor", "build",
                                                         "--hibf",
                                                         "--kmer 19",
                                                         "--window", std::to_string(window_size),
                                                         "--fpr 0.05",
                                                         "--threads", run_parallel ? "2" : "1",
                                                         "--output raptor.index",
                                                         pack_path(number_of_repeated_bins));
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);

    compare_index<raptor::index_structure::hibf>(ibf_path(number_of_repeated_bins,
                                                          window_size,
                                                          is_compressed::no,
                                                          is_hibf::yes),
                                                 "raptor.index");
}

TEST_P(build_hibf, with_shape)
{
    auto const [number_of_repeated_bins, window_size, run_parallel_tmp] = GetParam();
    bool const run_parallel = run_parallel_tmp && number_of_repeated_bins >= 32;

    cli_test_result const result = execute_app("raptor", "build",
                                                         "--hibf",
                                                         "--shape 1111111111111111111",
                                                         "--window", std::to_string(window_size),
                                                         "--fpr 0.05",
                                                         "--threads", run_parallel ? "2" : "1",
                                                         "--output raptor.index",
                                                         pack_path(number_of_repeated_bins));
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);

    compare_index<raptor::index_structure::hibf>(ibf_path(number_of_repeated_bins,
                                                          window_size,
                                                          is_compressed::no,
                                                          is_hibf::yes),
                                                 "raptor.index");
}

INSTANTIATE_TEST_SUITE_P(
    build_hibf_suite,
    build_hibf,
    testing::Combine(testing::Values(0, 16, 32), testing::Values(19, 23), testing::Values(true, false)),
    [] (testing::TestParamInfo<build_hibf::ParamType> const & info)
    {
        std::string name = std::to_string(std::max<int>(1, std::get<0>(info.param) * 4)) + "_bins_" +
                           std::to_string(std::get<1>(info.param)) + "_window_" +
                           (std::get<2>(info.param) ? "parallel" : "serial");
        return name;
    });

TEST_F(build_hibf, three_levels)
{
    cli_test_result const result = execute_app("raptor", "build",
                                                         "--hibf",
                                                         "--kmer 19",
                                                         "--window 19",
                                                         "--fpr 0.05",
                                                         "--threads 1",
                                                         "--output raptor.index",
                                                         data("three_levels.pack"));
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);

    compare_index<raptor::index_structure::hibf>(data("three_levels.hibf"), "raptor.index");
}

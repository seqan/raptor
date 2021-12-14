// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include "../cli_test.hpp"

struct search_hibf_preprocessing : public raptor_base,
                                   public testing::WithParamInterface<std::tuple<size_t, size_t, bool, size_t>> {};

TEST_P(search_hibf_preprocessing, pipeline)
{
    auto const [number_of_repeated_bins, window_size, run_parallel_tmp, number_of_errors] = GetParam();
    bool const run_parallel = run_parallel_tmp && number_of_repeated_bins >= 32;

    { // generate list of bins
        std::ofstream file{"raptor_cli_test.txt"};
        for (auto && file_path : get_repeated_bins(number_of_repeated_bins))
            file << file_path << '\n';
    }

    { // generate layout file
        std::ifstream input{pack_path(number_of_repeated_bins)};
        std::string line{};
        std::ofstream file{"raptor_cli_test.layout"};

        while (std::getline(input, line) && line.substr(0, 6) != "#FILES")
            file << line << '\n';
        file << line << '\n';
        while (std::getline(input, line))
        {
            std::string_view sv{line};
            size_t const filename_end = sv.find_first_of('\t');
            size_t const bin_name_end = filename_end - 3u;
            size_t const bin_name_start = sv.find_last_of('/', bin_name_end) + 1u;
            std::string_view bin_stem{sv.substr(bin_name_start, bin_name_end - bin_name_start)};
            file << "precomputed_minimisers/" << bin_stem << ".minimiser" << sv.substr(filename_end) << '\n';
        }
    }

    cli_test_result const result1 = execute_app("raptor", "build",
                                                          "--kmer 19",
                                                          "--window ", std::to_string(window_size),
                                                          "--compute-minimiser",
                                                          "--disable-cutoffs",
                                                          "--threads ", run_parallel ? "2" : "1",
                                                          "--output precomputed_minimisers",
                                                          "raptor_cli_test.txt");
    EXPECT_EQ(result1.out, std::string{});
    EXPECT_EQ(result1.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result1);

    cli_test_result const result2 = execute_app("raptor", "build",
                                                          "--hibf",
                                                          "--threads ", run_parallel ? "2" : "1",
                                                          "--output raptor.index",
                                                          "raptor_cli_test.layout");
    EXPECT_EQ(result2.out, std::string{});
    EXPECT_EQ(result2.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result2);

    compare_index<raptor::index_structure::hibf>(ibf_path(number_of_repeated_bins,
                                                          window_size,
                                                          is_compressed::no,
                                                          is_hibf::yes),
                                                 "raptor.index",
                                                 compare_extension::no);

    cli_test_result const result3 = execute_app("raptor", "search",
                                                          "--hibf",
                                                          "--output search.out",
                                                          "--error ", std::to_string(number_of_errors),
                                                          "--index ", "raptor.index",
                                                          "--query ", data("query.fq"));
    EXPECT_EQ(result3.out, std::string{});
    EXPECT_EQ(result3.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result3);

    compare_search(number_of_repeated_bins, number_of_errors, "search.out", is_empty::no, is_preprocessed::yes);
}

INSTANTIATE_TEST_SUITE_P(
    search_hibf_preprocessing_suite,
    search_hibf_preprocessing,
    testing::Combine(testing::Values(0, 16, 32),
                     testing::Values(19, 23),
                     testing::Values(true, false),
                     testing::Values(0, 1)),
    [] (testing::TestParamInfo<search_hibf_preprocessing::ParamType> const & info)
    {
        std::string name = std::to_string(std::max<int>(1, std::get<0>(info.param) * 4)) + "_bins_" +
                        std::to_string(std::get<1>(info.param)) + "_window_" +
                        (std::get<2>(info.param) ? "parallel" : "serial") +
                        std::to_string(std::get<3>(info.param)) + "_error";
        return name;
    });

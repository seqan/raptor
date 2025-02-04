// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <raptor/test/cli_test.hpp>
#include <raptor/update/dump_index.hpp>

struct build_hibf_layout : public raptor_base
{};

TEST_F(build_hibf_layout, pipeline)
{
    std::filesystem::path const data_filename = "raptor_cli_test.txt";
    std::filesystem::path const layout_filename = "raptor_cli_test.layout";
    std::string const index_filename = "raptor.index";
    std::filesystem::path const search_filename = "search.out";
#ifndef NDEBUG
    size_t const number_of_repeated_bins{2};
    size_t const initial_repeats{1};
#else
    size_t const number_of_repeated_bins{32};
    size_t const initial_repeats{16};
#endif
    size_t const number_of_errors{0}; // search

    { // generate sequence (data) input file
        std::ofstream file{data_filename};
        size_t dummy_group{}; // to avoid clustering by filenames in chopper sketch
        for (auto && file_path : get_repeated_bins(initial_repeats))
            file << file_path << '\t' << ++dummy_group << '\n';
    }

    ASSERT_TRUE(std::filesystem::exists(data_filename));

    { // build layout
        cli_test_result const result = execute_app("raptor",
                                                   "layout",
                                                   "--kmer 19",
                                                   "--threads 1",
                                                   "--input",
                                                   data_filename,
                                                   "--tmax 64",
                                                   "--empty-bin-fraction 0.001",
                                                   "--fpr 0.05",
                                                   "--output",
                                                   layout_filename);

        EXPECT_EQ(result.out, std::string{});
        EXPECT_EQ(result.err, std::string{});
        RAPTOR_ASSERT_ZERO_EXIT(result);
    }

    ASSERT_TRUE(std::filesystem::exists(layout_filename));

    std::cerr << std::filesystem::absolute(layout_filename) << '\n';

    { // build index
        cli_test_result const result = execute_app("raptor",
                                                   "build",
                                                   "--threads 1",
                                                   "--output",
                                                   index_filename + "_1",
                                                   "--quiet",
                                                   "--input",
                                                   layout_filename);

        EXPECT_EQ(result.out, std::string{});
        EXPECT_EQ(result.err, std::string{});
        RAPTOR_ASSERT_ZERO_EXIT(result);
    }

    size_t counter{1u};
    {
        for (auto && file_path : get_repeated_bins(number_of_repeated_bins - initial_repeats))
        {
            cli_test_result const result = execute_app("raptor",
                                                       "update",
                                                       "insert",
                                                       "--threads 1",
                                                       "--index",
                                                       index_filename + '_' + std::to_string(counter),
                                                       "--output",
                                                       index_filename + '_' + std::to_string(counter + 1u),
                                                       "--insert",
                                                       file_path);

            // EXPECT_EQ(result.out, std::string{});
            EXPECT_EQ(result.err, std::string{});
            RAPTOR_ASSERT_ZERO_EXIT(result);
            ++counter;
        }
    }

    {
        std::ifstream index_file{index_filename + '_' + std::to_string(counter)};
        cereal::BinaryInputArchive archive{index_file};
        raptor::raptor_index<raptor::index_structure::hibf> index;
        archive(index);
        raptor::dump_index(index);
        // for (auto & ibf : index.ibf().ibf_vector)
        // {
        //     std::cerr << ibf.bin_count() << '\n';
        //     for (auto val : ibf.occupancy)
        //     {
        //         std::cout << val << ',';
        //     }
        //     std::cout << '\n';
        // }
    }

    { // check with search if index contains expected input
        cli_test_result const result = execute_app("raptor",
                                                   "search",
                                                   "--output",
                                                   search_filename,
                                                   "--error",
                                                   std::to_string(number_of_errors),
                                                   "--index",
                                                   index_filename + '_' + std::to_string(counter),
                                                   "--quiet",
                                                   "--query",
                                                   data("query.fq"));
        EXPECT_EQ(result.out, std::string{});
        EXPECT_EQ(result.err, std::string{});
        RAPTOR_ASSERT_ZERO_EXIT(result);
    }

    compare_search(number_of_repeated_bins, number_of_errors, search_filename.c_str());
}

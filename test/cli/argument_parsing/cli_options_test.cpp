// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <seqan3/test/snippet/create_temporary_snippet_file.hpp>

#include "../cli_test.hpp"

struct argparse_build : public raptor_base {};
struct argparse_main : public raptor_base {};
struct argparse_search : public raptor_base {};
struct argparse_upgrade : public raptor_base {};

seqan3::test::create_temporary_snippet_file tmp_index_file{"tmp.index", "\nsome_content"};
seqan3::test::create_temporary_snippet_file dummy_sequence_file{"dummy.fasta", "\n>ID\nACGTC"};
seqan3::test::create_temporary_snippet_file tmp_bin_list_file{"all_bins.txt", std::string{"\n"} +
                                                                              dummy_sequence_file.file_path.string()};

TEST_F(argparse_main, no_options)
{
    cli_test_result const result = execute_app("raptor");
    std::string const expected
    {
        "Raptor - A fast and space-efficient pre-filter for querying very large collections of nucleotide sequences.\n"
        "===========================================================================================================\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);
}

TEST_F(argparse_build, no_options)
{
    cli_test_result const result = execute_app("raptor", "build");
    std::string const expected
    {
        "Raptor - A fast and space-efficient pre-filter for querying very large collections of nucleotide sequences.\n"
        "===========================================================================================================\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);
}

TEST_F(argparse_search, no_options)
{
    cli_test_result const result = execute_app("raptor", "search");
    std::string const expected
    {
        "Raptor - A fast and space-efficient pre-filter for querying very large collections of nucleotide sequences.\n"
        "===========================================================================================================\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);
}

TEST_F(argparse_main, no_subparser)
{
    cli_test_result const result = execute_app("raptor", "foo");
    std::string const expected
    {
        "[Error] You either forgot or misspelled the subcommand! Please specify which sub-program you want to use: one "
        "of [build,search,socks,upgrade]. Use -h/--help for more information.\n"
    };
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_main, unknown_option)
{
    cli_test_result const result = execute_app("raptor", "-v");
    std::string const expected
    {
        "[Error] You either forgot or misspelled the subcommand! Please specify which sub-program you want to use: one "
        "of [build,search,socks,upgrade]. Use -h/--help for more information.\n"
    };
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, input_missing)
{
    cli_test_result const result = execute_app("raptor", "build",
                                                         "--size 8m",
                                                         "--output ./index.raptor");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Not enough positional arguments provided (Need at least 1). See "
                                      "-h/--help for more information.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, input_invalid)
{
    cli_test_result const result = execute_app("raptor", "build",
                                                         "--size 8m",
                                                         "--output ./index.raptor",
                                                         "nonexistent");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Validation failed for positional option 1: The file \"nonexistent\" does"
                                      " not exist!\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, output_missing)
{
    cli_test_result const result = execute_app("raptor", "build",
                                                         "--size 8m",
                                                         tmp_bin_list_file.file_path);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Option --output is required but not set.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, directory_missing)
{
    cli_test_result const result = execute_app("raptor", "build",
                                                         "--size 8m",
                                                         "--compute-minimiser",
                                                         tmp_bin_list_file.file_path);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Option --output is required but not set.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, alias)
{
    cli_test_result const result = execute_app("raptor", "build",
                                                         "--size 8m",
                                                         "--compute-minimizer",
                                                         tmp_bin_list_file.file_path);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Option --output is required but not set.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, size_missing)
{
    cli_test_result const result = execute_app("raptor", "build",
                                                         "--output ./index.raptor",
                                                         tmp_bin_list_file.file_path);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Option --size is required but not set.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, size_wrong_space)
{
    cli_test_result const result = execute_app("raptor", "build",
                                                         "--size 8 m",
                                                         "--output ./index.raptor",
                                                         tmp_bin_list_file.file_path);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Validation failed for option --size: Value 8 must be an integer "
                                      "followed by [k,m,g,t] (case insensitive).\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, size_wrong_suffix)
{
    cli_test_result const result = execute_app("raptor", "build",
                                                         "--size 8x",
                                                         "--output index.raptor",
                                                         tmp_bin_list_file.file_path);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Validation failed for option --size: Value 8x must be an integer "
                                      "followed by [k,m,g,t] (case insensitive).\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, kmer_window)
{
    cli_test_result const result = execute_app("raptor", "build",
                                                         "--kmer 20",
                                                         "--window 19",
                                                         "--size 8m",
                                                         "--output index.raptor",
                                                         tmp_bin_list_file.file_path);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] The k-mer size cannot be bigger than the window size.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, kmer_shape)
{
    cli_test_result const result = execute_app("raptor", "build",
                                                         "--kmer 20",
                                                         "--shape 11",
                                                         "--size 8m",
                                                         "--output index.raptor",
                                                         tmp_bin_list_file.file_path);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] You cannot set both shape and k-mer arguments.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, zero_threads)
{
    cli_test_result const result = execute_app("raptor", "build",
                                                         "--kmer 20",
                                                         "--threads 0",
                                                         "--size 8m",
                                                         "--output index.raptor",
                                                         tmp_bin_list_file.file_path);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Validation failed for option --threads: The value must be a positive "
                                      "integer.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, no_bins_in_file)
{
    seqan3::test::create_temporary_snippet_file tmp_bin_list_empty{"empty.txt", std::string{"\n"}};

    cli_test_result const result = execute_app("raptor", "build",
                                                         "--kmer 20",
                                                         "--size 8m",
                                                         "--output index.raptor",
                                                         tmp_bin_list_empty.file_path);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] The list of input files cannot be empty.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, empty_file_in_bin)
{
    seqan3::test::create_temporary_snippet_file empty_sequence_file{"empty.fasta", std::string{"\n"}};
    seqan3::test::create_temporary_snippet_file tmp_empty_bin_file{"empty_bin.txt",
                                                                   std::string{"\n"} +
                                                                   empty_sequence_file.file_path.string()};

    cli_test_result const result = execute_app("raptor", "build",
                                                         "--kmer 20",
                                                         "--size 8m",
                                                         "--output index.raptor",
                                                         tmp_empty_bin_file.file_path);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] The file " + tmp_empty_bin_file.file_path.parent_path().string() +
                                      "/empty.fasta is empty.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, mixed_input)
{
    seqan3::test::create_temporary_snippet_file minimiser_file{"bin1.minimiser", std::string{"\n0"}};
    seqan3::test::create_temporary_snippet_file mixed_bin_file{"mixed.txt", std::string{"\n"} +
                                                                            minimiser_file.file_path.string() +
                                                                            std::string{"\nbin2.fasta"}};

    cli_test_result const result = execute_app("raptor", "build",
                                                         "--kmer 20",
                                                         "--size 8m",
                                                         "--output index.raptor",
                                                         mixed_bin_file.file_path);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] You cannot mix sequence and minimiser files as input.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, wrong_parts)
{
    cli_test_result const result = execute_app("raptor", "build",
                                                         "--kmer 20",
                                                         "--size 8m",
                                                         "--parts 3",
                                                         "--output index.raptor",
                                                         tmp_bin_list_file.file_path);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Validation failed for option --parts: The value must be a power of "
                                      "two.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_search, ibf_missing)
{
    cli_test_result const result = execute_app("raptor", "search",
                                                         "--query ", data("query.fq"),
                                                         "--output search.out");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Option --index is required but not set.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_search, ibf_wrong)
{
    cli_test_result const result = execute_app("raptor", "search",
                                                         "--query ", data("query.fq"),
                                                         "--index foo.index",
                                                         "--output search.out");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] The file \"foo.index\" does not exist!\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_search, query_missing)
{
    cli_test_result const result = execute_app("raptor", "search",
                                                         "--index ", tmp_index_file.file_path,
                                                         "--output search.out");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Option --query is required but not set.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_search, query_wrong)
{
    cli_test_result const result = execute_app("raptor", "search",
                                                         "--query foo.fasta",
                                                         "--index ", tmp_index_file.file_path,
                                                         "--output search.out");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Validation failed for option --query: The file \"foo.fasta\" does not "
                                      "exist!\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_search, output_missing)
{
    cli_test_result const result = execute_app("raptor", "search",
                                                         "--query ", data("query.fq"),
                                                         "--index ", tmp_index_file.file_path);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Option --output is required but not set.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_search, old_index)
{
    cli_test_result const result = execute_app("raptor", "search",
                                                         "--query ", data("query.fq"),
                                                         "--index ", data("1_1.index"),
                                                         "--output search.out");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Unsupported index version. Check raptor upgrade.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_upgrade, kmer_window)
{
    cli_test_result const result = execute_app("raptor", "upgrade",
                                                         "--bins ", tmp_bin_list_file.file_path,
                                                         "--input ", data("1_1.index"),
                                                         "--output index.raptor",
                                                         "--window 19",
                                                         "--kmer 20");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] The k-mer size cannot be bigger than the window size.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <seqan3/test/tmp_directory.hpp>

#include <raptor/test/cli_test.hpp>

struct argparse_build : public raptor_base
{};
struct argparse_main : public raptor_base
{};
struct argparse_search : public raptor_base
{};
struct argparse_upgrade : public raptor_base
{};

seqan3::test::tmp_directory const test_directory{};
std::filesystem::path const tmp_index_file = []()
{
    std::filesystem::path const tmp_index_file = test_directory.path() / "tmp.index";
    std::ofstream file{tmp_index_file};
    file << "some_content";
    return tmp_index_file;
}();
std::filesystem::path const dummy_sequence_file = []()
{
    std::filesystem::path const dummy_sequence_file = test_directory.path() / "dummy.fasta";
    std::ofstream file{dummy_sequence_file};
    file << ">ID\nACGTCACGATCGTACGATCGATCGATCG";
    return dummy_sequence_file;
}();
std::filesystem::path const tmp_bin_list_file = [](std::filesystem::path const dummy_sequence_file)
{
    std::filesystem::path const tmp_bin_list_file = test_directory.path() / "all_bins.txt";
    std::ofstream file{tmp_bin_list_file};
    file << dummy_sequence_file.c_str();
    return tmp_bin_list_file;
}(dummy_sequence_file);

std::filesystem::path const tmp_bin_list_empty = []()
{
    std::filesystem::path const tmp_bin_list_empty = test_directory.path() / "empty.txt";
    std::ofstream file{tmp_bin_list_empty};
    return tmp_bin_list_empty;
}();

std::filesystem::path const empty_sequence_file = []()
{
    std::filesystem::path const empty_sequence_file = test_directory.path() / "empty.fasta";
    std::ofstream file{empty_sequence_file};
    return empty_sequence_file;
}();

std::filesystem::path const tmp_empty_bin_file = [](std::filesystem::path const empty_sequence_file)
{
    std::filesystem::path const tmp_empty_bin_file = test_directory.path() / "empty_bin.txt";
    std::ofstream file{tmp_empty_bin_file};
    file << empty_sequence_file.c_str();
    return tmp_empty_bin_file;
}(empty_sequence_file);

std::filesystem::path const minimiser_file = []()
{
    std::filesystem::path const minimiser_file = test_directory.path() / "bin1.minimiser";
    std::ofstream file{minimiser_file};
    file << "0";
    return minimiser_file;
}();

std::filesystem::path const mixed_bin_file = [](std::filesystem::path const minimiser_file)
{
    std::filesystem::path const mixed_bin_file = test_directory.path() / "mixed.txt";
    std::ofstream file{mixed_bin_file};
    file << minimiser_file.c_str() << "\nbin2.fasta";
    return mixed_bin_file;
}(minimiser_file);

TEST_F(argparse_main, no_options)
{
    cli_test_result const result = execute_app("raptor");
    std::string const expected{
        "Raptor - A fast and space-efficient pre-filter for querying very large collections of nucleotide sequences.\n"
        "===========================================================================================================\n"
        "    Try -h or --help for more information.\n"};
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);
}

TEST_F(argparse_build, no_options)
{
    cli_test_result const result = execute_app("raptor", "build");
    std::string const expected{"Raptor-build - A fast and space-efficient pre-filter for querying very large "
                               "collections of nucleotide sequences.\n"
                               "======================================================================================="
                               "==========================\n"
                               "    Try -h or --help for more information.\n"};
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);
}

TEST_F(argparse_search, no_options)
{
    cli_test_result const result = execute_app("raptor", "search");
    std::string const expected{"Raptor-search - A fast and space-efficient pre-filter for querying very large "
                               "collections of nucleotide sequences.\n"
                               "======================================================================================="
                               "===========================\n"
                               "    Try -h or --help for more information.\n"};
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);
}

TEST_F(argparse_main, no_subparser)
{
    cli_test_result const result = execute_app("raptor", "foo");
    std::string const expected{
        "[Error] You either forgot or misspelled the subcommand! Please specify which sub-program you want to use: one "
        "of [build, layout, prepare, search, upgrade]. Use -h/--help for more information.\n"};
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_main, unknown_option)
{
    cli_test_result const result = execute_app("raptor", "-v");
    std::string const expected{
        "[Error] You either forgot or misspelled the subcommand! Please specify which sub-program you want to use: one "
        "of [build, layout, prepare, search, upgrade]. Use -h/--help for more information.\n"};
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, input_missing)
{
    cli_test_result const result = execute_app("raptor", "build", "--output ./index.raptor");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err,
              std::string{"[Error] Not enough positional arguments provided (Need at least 1). See "
                          "-h/--help for more information.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, input_invalid)
{
    cli_test_result const result = execute_app("raptor", "build", "--output ./index.raptor", "nonexistent");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err,
              std::string{"[Error] Validation failed for positional option 1: The file \"nonexistent\" does"
                          " not exist!\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, output_missing)
{
    cli_test_result const result = execute_app("raptor", "build", tmp_bin_list_file);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Option --output is required but not set.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, directory_missing)
{
    cli_test_result const result = execute_app("raptor", "build", "--compute-minimiser", tmp_bin_list_file);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Option --output is required but not set.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, alias)
{
    cli_test_result const result = execute_app("raptor", "build", "--compute-minimizer", tmp_bin_list_file);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Option --output is required but not set.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, kmer_window)
{
    cli_test_result const result =
        execute_app("raptor", "build", "--kmer 20", "--window 19", "--output index.raptor", tmp_bin_list_file);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] The k-mer size cannot be bigger than the window size.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, kmer_shape)
{
    cli_test_result const result =
        execute_app("raptor", "build", "--kmer 20", "--shape 11", "--output index.raptor", tmp_bin_list_file);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] You cannot set both shape and k-mer arguments.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, zero_threads)
{
    cli_test_result const result =
        execute_app("raptor", "build", "--kmer 20", "--threads 0", "--output index.raptor", tmp_bin_list_file);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err,
              std::string{"[Error] Validation failed for option --threads: The value must be a positive "
                          "integer.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, no_bins_in_file)
{
    cli_test_result const result =
        execute_app("raptor", "build", "--kmer 20", "--output index.raptor", tmp_bin_list_empty);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] The list of input files cannot be empty.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, empty_file_in_bin)
{
    cli_test_result const result =
        execute_app("raptor", "build", "--kmer 20", "--output index.raptor", tmp_empty_bin_file);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(
        result.err,
        std::string{"[Error] The file " + tmp_empty_bin_file.parent_path().string() + "/empty.fasta is empty.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, mixed_input)
{
    cli_test_result const result = execute_app("raptor", "build", "--kmer 20", "--output index.raptor", mixed_bin_file);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] You cannot mix sequence and minimiser files as input.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, wrong_parts)
{
    cli_test_result const result =
        execute_app("raptor", "build", "--kmer 20", "--parts 3", "--output index.raptor", tmp_bin_list_file);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err,
              std::string{"[Error] Validation failed for option --parts: The value must be a power of "
                          "two.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_search, ibf_missing)
{
    cli_test_result const result = execute_app("raptor", "search", "--query ", data("query.fq"), "--output search.out");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Option --index is required but not set.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_search, ibf_wrong)
{
    cli_test_result const result =
        execute_app("raptor", "search", "--query ", data("query.fq"), "--index foo.index", "--output search.out");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] The file \"foo.index\" does not exist!\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_search, query_missing)
{
    cli_test_result const result = execute_app("raptor", "search", "--index ", tmp_index_file, "--output search.out");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Option --query is required but not set.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_search, query_wrong)
{
    cli_test_result const result =
        execute_app("raptor", "search", "--query foo.fasta", "--index ", tmp_index_file, "--output search.out");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err,
              std::string{"[Error] Validation failed for option --query: The file \"foo.fasta\" does not "
                          "exist!\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_search, output_missing)
{
    cli_test_result const result =
        execute_app("raptor", "search", "--query ", data("query.fq"), "--index ", tmp_index_file);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Option --output is required but not set.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_search, old_index)
{
    cli_test_result const result = execute_app("raptor",
                                               "search",
                                               "--query ",
                                               data("query.fq"),
                                               "--index ",
                                               data("1_1.index"),
                                               "--output search.out");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Unsupported index version. Check raptor upgrade.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_search, error_treshold)
{
    cli_test_result const result = execute_app("raptor",
                                               "search",
                                               "--query ",
                                               data("query.fq"),
                                               "--index ",
                                               data("1bins23window.index"),
                                               "--error 0",
                                               "--threshold 0.4",
                                               "--output search.out");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] You cannot set both error and threshold arguments.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_search, empty_query)
{
    cli_test_result const result = execute_app("raptor",
                                               "search",
                                               "--query ",
                                               data("empty.fq"),
                                               "--index ",
                                               data("1bins23window.index"),
                                               "--output search.out");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] The query file is empty.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_search, queries_too_short)
{
    cli_test_result const result = execute_app("raptor",
                                               "search",
                                               "--query ",
                                               data("too_short.fq"),
                                               "--index ",
                                               data("1bins23window.index"),
                                               "--output search.out");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] The query size (21) is too short to be used with window size 23.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_upgrade, kmer_window)
{
    cli_test_result const result = execute_app("raptor",
                                               "upgrade",
                                               "--bins ",
                                               tmp_bin_list_file,
                                               "--input ",
                                               data("1_1.index"),
                                               "--output index.raptor",
                                               "--window 19",
                                               "--kmer 20");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] The k-mer size cannot be bigger than the window size.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

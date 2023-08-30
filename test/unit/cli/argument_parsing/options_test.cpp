// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <chopper/prefixes.hpp>

#include <raptor/test/cli_test.hpp>
#include <raptor/test/tmp_test_file.hpp>

struct argparse_build : public raptor_base
{};
struct argparse_main : public raptor_base
{};
struct argparse_search : public raptor_base
{};
struct argparse_upgrade : public raptor_base
{};
struct argparse_prepare : public raptor_base
{};
struct argparse_layout : public raptor_base
{};

raptor::test::tmp_test_file const test_files{};
// clang-format off
std::filesystem::path const
    tmp_index_file         = test_files.create("tmp.index",     "some_content"),
    dummy_sequence_file    = test_files.create("dummy.fasta",   ">ID\nACGTCACGATCGTACGATCGATCGATCG"),
    tmp_bin_list_file      = test_files.create("all_bins.txt",   dummy_sequence_file.c_str()),
    tmp_bin_list_empty     = test_files.create("empty.txt"),
    tmp_bin_list_no_exist  = test_files.create("no_exist.txt", "raptor_options_test_non_existing_file.fa"),
    tmp_bin_list_corrupted = test_files.create("corrupted.txt",  chopper::prefix::meta_chopper_user_bins_start, '\n', chopper::prefix::meta_chopper_user_bins_end),
    empty_sequence_file    = test_files.create("empty.fasta"),
    tmp_empty_bin_file     = test_files.create("empty_bin.txt",  empty_sequence_file.c_str()),
    header_file            = test_files.create("bin1.header",    "1111111111111111111\t19\t0\t1\n"),
    minimiser_file         = test_files.create("bin1.minimiser", "0"),
    minimiser_list         = test_files.create("minimiser.list", minimiser_file.c_str()),
    mixed_bin_file         = test_files.create("mixed.txt",      minimiser_file.c_str(), "\nbin2.fasta"),
    old_layout_file        = test_files.create("old.layout",     "#TOP_LEVEL_IBF");
// clang-format on

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
    std::string const expected{
        "Raptor-build - A fast and space-efficient pre-filter for querying very large collections of nucleotide sequen"
        "ces.\n======================================================================================================="
        "==========\n    raptor build --input <file> --output <file> [--threads <number>] [--quiet]\n    [--kmer <numb"
        "er>|--shape <01-pattern>] [--window <number>] [--fpr\n    <number>] [--hash <number>] [--parts <number>]\n   "
        " Try -h or --help for more information.\n"};
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);
}

TEST_F(argparse_search, no_options)
{
    cli_test_result const result = execute_app("raptor", "search");
    std::string const expected{
        "Raptor-search - A fast and space-efficient pre-filter for querying very large collections of nucleotide seque"
        "nces.\n======================================================================================================"
        "============\n    raptor search --index <file> --query <file> --output <file> [--threads\n    <number>] [--qu"
        "iet] [--error <number>|--threshold <number>]\n    [--query_length <number>] [--tau <number>] [--pmax <number>"
        "]\n    [--cache-thresholds]\n    Try -h or --help for more information.\n"};
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);
}

TEST_F(argparse_main, no_subparser)
{
    cli_test_result const result = execute_app("raptor", "foo");
    std::string const expected{
        "[Error] You misspelled the subcommand! Please specify which sub-program you want to use: one "
        "of [build, layout, prepare, search, upgrade]. Use -h/--help for more information.\n"};
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_main, unknown_option)
{
    cli_test_result const result = execute_app("raptor", "-v");
    std::string const expected{
        "[Error] You misspelled the subcommand! Please specify which sub-program you want to use: one "
        "of [build, layout, prepare, search, upgrade]. Use -h/--help for more information.\n"};
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_layout, no_bins_in_file)
{
    cli_test_result const result = execute_app("raptor", "layout", "--input-file", tmp_bin_list_empty);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] The input file is empty.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, input_missing)
{
    cli_test_result const result = execute_app("raptor", "build", "--output ./index.raptor");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Option --input is required but not set.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, input_invalid)
{
    cli_test_result const result = execute_app("raptor", "build", "--output ./index.raptor", "--input nonexistent");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err,
              std::string{"[Error] Validation failed for option --input: The file \"nonexistent\" does not exist!\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, output_missing)
{
    cli_test_result const result = execute_app("raptor", "build", "--input", tmp_bin_list_file);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Option --output is required but not set.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, kmer_window)
{
    cli_test_result const result = execute_app("raptor",
                                               "build",
                                               "--kmer 20",
                                               "--window 19",
                                               "--output index.raptor",
                                               "--input",
                                               tmp_bin_list_file);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] The k-mer size cannot be bigger than the window size.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, kmer_shape)
{
    cli_test_result const result = execute_app("raptor",
                                               "build",
                                               "--kmer 20",
                                               "--shape 11",
                                               "--output index.raptor",
                                               "--input",
                                               tmp_bin_list_file);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] You cannot set both shape and k-mer arguments.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, zero_threads)
{
    cli_test_result const result = execute_app("raptor",
                                               "build",
                                               "--kmer 20",
                                               "--threads 0",
                                               "--output index.raptor",
                                               "--input",
                                               tmp_bin_list_file);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err,
              std::string{"[Error] Validation failed for option --threads: The value must be a positive "
                          "integer.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, no_bins_in_file)
{
    cli_test_result const result =
        execute_app("raptor", "build", "--kmer 20", "--output index.raptor", "--input", tmp_bin_list_empty);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] The input file is empty.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, bin_does_not_exist)
{
    cli_test_result const result =
        execute_app("raptor", "build", "--kmer 20", "--output index.raptor", "--input", tmp_bin_list_no_exist);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] The file raptor_options_test_non_existing_file.fa does not exist.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, corrupted_input)
{
    cli_test_result const result =
        execute_app("raptor", "build", "--kmer 20", "--output index.raptor", "--input", tmp_bin_list_corrupted);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] The list of input files cannot be empty.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, empty_file_in_bin)
{
    cli_test_result const result =
        execute_app("raptor", "build", "--kmer 20", "--output index.raptor", "--input", tmp_empty_bin_file);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(
        result.err,
        std::string{"[Error] The file " + tmp_empty_bin_file.parent_path().string() + "/empty.fasta is empty.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, mixed_input)
{
    cli_test_result const result =
        execute_app("raptor", "build", "--kmer 20", "--output index.raptor", "--input", mixed_bin_file);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] You cannot mix sequence and minimiser files as input.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, wrong_parts)
{
    cli_test_result const result =
        execute_app("raptor", "build", "--kmer 20", "--parts 3", "--output index.raptor", "--input", tmp_bin_list_file);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err,
              std::string{"[Error] Validation failed for option --parts: The value must be a power of "
                          "two.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, partitioned_parts)
{
    cli_test_result const result =
        execute_app("raptor", "build", "--parts 4", "--output index.raptor", "--input", data("three_levels.pack"));
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] The HIBF cannot yet be partitioned.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, minimiser_and_shape)
{
    cli_test_result const result =
        execute_app("raptor", "build", "--shape 11111", "--output index.raptor", "--input", minimiser_list);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] You cannot set --shape when using minimiser files as input.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, minimiser_and_kmer)
{
    cli_test_result const result =
        execute_app("raptor", "build", "--kmer 19", "--output index.raptor", "--input", minimiser_list);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] You cannot set --kmer when using minimiser files as input.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, minimiser_and_window)
{
    cli_test_result const result =
        execute_app("raptor", "build", "--window 19", "--output index.raptor", "--input", minimiser_list);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] You cannot set --window when using minimiser files as input.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, layout_config_and_options)
{
    cli_test_result const result = execute_app("raptor",
                                               "build",
                                               "--kmer 20",
                                               "--hash 3",
                                               "--fpr 0.01",
                                               "--output index.raptor",
                                               "--quiet",
                                               "--input",
                                               data("test.layout"));
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(
        result.err,
        std::string{"[WARNING] Given k-mer size(20) differs from k-mer size in the layout file (19). The results may "
                    "be suboptimal. If this was a conscious decision, you can ignore this warning.\n[WARNING] Given "
                    "hash function count (3) differs from hash function count in the layout file (2). The results may "
                    "be suboptimal. If this was a conscious decision, you can ignore this warning.\n[WARNING] Given "
                    "false positive rate (0.01) differs from false positive rate in the layout file (0.05). The "
                    "results may be suboptimal. If this was a conscious decision, you can ignore this warning."});
    RAPTOR_ASSERT_ZERO_EXIT(result);
}

TEST_F(argparse_build, layout_too_small_window)
{
    cli_test_result const result =
        execute_app("raptor", "build", "--window 18", "--output index.raptor", "--input", data("test.layout"));
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] The k-mer size cannot be bigger than the window size.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_build, layout_old)
{
    cli_test_result const result = execute_app("raptor", "build", "--output index.raptor", "--input", old_layout_file);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err,
              std::string{"[Error] The input file was determined to be a layout. However, the first line starts "
                          "with '#' (old format) instead of '@' (new format).\n"});
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
                                               data("2.0.index"),
                                               "--output search.out");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Unsupported index version. Check raptor upgrade.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_search, ambiguous_index)
{
    {
        std::ofstream monolithic_index{"raptor.index"};
        std::ofstream partitioned_index{"raptor.index_0"};
    }

    cli_test_result const result =
        execute_app("raptor", "search", "--query ", data("query.fq"), "--index raptor.index", "--output search.out");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err,
              std::string{"[Error] Ambiguous index. Both monolithic (raptor.index) and partitioned index "
                          "(raptor.index_0) exist. Please rename the monolithic index.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_search, compressed_index)
{

    cli_test_result const result = execute_app("raptor",
                                               "search",
                                               "--query ",
                                               data("query.fq"),
                                               "--index",
                                               cli_test::data("64bins23windowc.index"),
                                               "--output search.out");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Cannot read index: Index cannot be compressed.\n"});
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
    EXPECT_EQ(result.err,
              std::string{"[Error] The (minimal) query length (21) is too short to be used with window size 23.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_search, queries_have_variance)
{
    cli_test_result const result = execute_app("raptor",
                                               "search",
                                               "--query ",
                                               data("query_variance.fq"),
                                               "--index ",
                                               data("1bins23window.index"),
                                               "--quiet",
                                               "--output search.out");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(
        result.err,
        std::string{
            "[WARNING] There is variance in the provided queries. The shortest length is 28. The longest length is 65. "
            "The tresholding will use a single query length (65). Therefore, results may be inprecise.\n"});
    RAPTOR_ASSERT_ZERO_EXIT(result);
}

TEST_F(argparse_search, queries_unsupported_length)
{
    std::filesystem::path const query_file = test_files.path() / "unsupported_length.fa";
    {
        std::ofstream os{query_file};
        // Add a long query. 65536 = 2^16, which is exactly one more than uint16_t can represent.
        os << ">query1\n";
        os << std::string(65536, 'A') << '\n';
        // Add a query shorter than the window size such that we throw.
        // Searching the long query takes too long for testing.
        os << ">query2\n";
        os << std::string(22, 'A') << '\n';
    }

    cli_test_result const result = execute_app("raptor",
                                               "search",
                                               "--query ",
                                               query_file,
                                               "--index ",
                                               data("1bins23window.index"),
                                               "--output search.out");

    std::string cerr_message{};
    cerr_message.reserve(553);

    // First part: We have a short and a long query -> high variance.
    cerr_message +=
        "[WARNING] There is variance in the provided queries. The shortest length is 22. The longest length is 65536. "
        "The tresholding will use a single query length (65536). Therefore, results may be inprecise.\n";
    // Second part: The actual warning for our query being too long.
    cerr_message += "[WARNING] There are queries which exceed the maximum safely supported length of 65535. Results "
                    "may be wrong, especially when using window size == k-mer size. If you need longer queries to be "
                    "supported, please open an issue at https://github.com/seqan/raptor/issues.\n";
    // Third part: Error message for the short query being too short. Early exit because querying takes too long.
    cerr_message += "[Error] The (minimal) query length (22) is too short to be used with window size 23.\n";
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, cerr_message);
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_upgrade, exclusive_options)
{
    {
        std::ofstream monolithic_index{"raptor.index"};
        std::ofstream bins{"bins.list"};
    }

    cli_test_result const result = execute_app("raptor",
                                               "upgrade",
                                               "--input raptor.index",
                                               "--output upgrade.index",
                                               "--fpr 0.05",
                                               "--bins bins.list");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] You cannot set both --fpr and --bins.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_upgrade, ambiguous_index)
{
    {
        std::ofstream monolithic_index{"raptor.index"};
        std::ofstream partitioned_index{"raptor.index_0"};
    }

    cli_test_result const result =
        execute_app("raptor", "upgrade", "--input raptor.index", "--output upgrade.index", "--fpr 0.05");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err,
              std::string{"[Error] Ambiguous index. Both monolithic (raptor.index) and partitioned index "
                          "(raptor.index_0) exist. Please rename the monolithic index.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_upgrade, unsupported_index)
{
    cli_test_result const result = execute_app("raptor",
                                               "upgrade",
                                               "--input",
                                               data("1bins23window.index"),
                                               "--output upgrade.index",
                                               "--fpr 0.05");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Unsupported index version. Use Raptor 2.0's upgrade first.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(argparse_prepare, no_options)
{
    cli_test_result const result = execute_app("raptor", "prepare");
    std::string const expected{
        "Raptor-prepare - A fast and space-efficient pre-filter for querying very large collections of nucleotide sequ"
        "ences.\n====================================================================================================="
        "==============\n    raptor prepare --input <file> --output <directory> [--threads <number>]\n    [--quiet] [-"
        "-kmer <number>|--shape <01-pattern>] [--window <number>]\n    [--kmer-count-cutoff <number>|--use-filesize-de"
        "pendent-cutoff]\n    Try -h or --help for more information.\n"};
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);
}

TEST_F(argparse_prepare, cutoffs)
{
    cli_test_result const result = execute_app("raptor",
                                               "prepare",
                                               "--output directory",
                                               "--kmer-count-cutoff 1",
                                               "--use-filesize-dependent-cutoff",
                                               "--input",
                                               tmp_bin_list_file);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err,
              std::string{"[Error] You cannot use both --kmer-count-cutoff and --use-filesize-dependent-cutoff.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

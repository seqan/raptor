#include <seqan3/std/ranges>     // range comparisons
#include <string>                // strings
#include <vector>                // vectors

#include "cli_test.hpp"

TEST_F(raptor, no_options)
{
    cli_test_result result = execute_app("raptor");
    std::string expected
    {
        "Raptor - A fast and space-efficient pre-filter for querying very large collections of nucleotide sequences.\n"
        "===========================================================================================================\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(raptor_build, no_options)
{
    cli_test_result result = execute_app("raptor", "build");
    std::string expected
    {
        "Raptor - A fast and space-efficient pre-filter for querying very large collections of nucleotide sequences.\n"
        "===========================================================================================================\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(raptor_search, no_options)
{
    cli_test_result result = execute_app("raptor", "search");
    std::string expected
    {
        "Raptor - A fast and space-efficient pre-filter for querying very large collections of nucleotide sequences.\n"
        "===========================================================================================================\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(raptor, no_subparser)
{
    cli_test_result result = execute_app("raptor", "foo");
    std::string expected
    {
        "[Error] Please specify which sub program you want to use (one of [build,search]). Use -h/--help for more information.\n"
    };
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(raptor, unknown_option)
{
    cli_test_result result = execute_app("raptor", "-v");
    std::string expected
    {
        "[Error] Please specify which sub program you want to use (one of [build,search]). Use -h/--help for more information.\n"
    };
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(raptor_build, output_missing)
{
    cli_test_result result = execute_app("raptor", "build", "--size 8m", expand_bins(64));
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Option --output is required but not set.\n"});
}

TEST_F(raptor_build, output_wrong)
{
    cli_test_result result = execute_app("raptor", "build", "--size 8m", "--output foo/out.ibf", expand_bins(64));
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Cannot write \"foo/out.ibf\"!\n"});
}

TEST_F(raptor_build, directory_missing)
{
    cli_test_result result = execute_app("raptor", "build", "--size 8m", "--compute-minimiser", expand_bins(64));
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Option --output is required but not set.\n"});
}

TEST_F(raptor_build, directory_wrong)
{
    cli_test_result result = execute_app("raptor", "build", "--size 8m", "--compute-minimiser", "--output foo/bar", expand_bins(64));
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Cannot create directory: \"foo/bar\"!\n"});
}

TEST_F(raptor_build, size_missing)
{
    cli_test_result result = execute_app("raptor", "build", "--output ./ibf.out", expand_bins(64));
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Option --size is required but not set.\n"});
}

TEST_F(raptor_build, size_wrong_space)
{
    cli_test_result result = execute_app("raptor", "build", "--size 8 m", "--output ./ibf.out", expand_bins(64));
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Validation failed for option --size: Value 8 must be an integer followed by [k,m,g,t] (case insensitive).\n"});
}

TEST_F(raptor_build, size_wrong_suffix)
{
    cli_test_result result = execute_app("raptor", "build", "--size 8x", "--output ibf.out", expand_bins(64));
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Validation failed for option --size: Value 8x must be an integer followed by [k,m,g,t] (case insensitive).\n"});
}

TEST_F(raptor_build, kmer_window)
{
    cli_test_result result = execute_app("raptor", "build", "--kmer 20", "--window 19", "--size 8m", "--output ibf.out", expand_bins(64));
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] The k-mer size cannot be bigger than the window size.\n"});
}

TEST_F(raptor_search, ibf_missing)
{
    cli_test_result result = execute_app("raptor", "search", "--query ", data("example_data/64/reads/all.fastq"), "--output search.out");
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Option --index is required but not set.\n"});
}

TEST_F(raptor_search, ibf_wrong)
{
    cli_test_result result = execute_app("raptor", "search", "--query ", data("example_data/64/reads/all.fastq"), "--index foo.ibf", "--output search.out");
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Validation failed for option --index: The file \"foo.ibf\" does not exist!\n"});
}

TEST_F(raptor_search, query_missing)
{
    cli_test_result result = execute_app("raptor", "search", "--index ", data("expected_results/b64_k19_w23_s8m.ibf"), "--output search.out");
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Option --query is required but not set.\n"});
}

TEST_F(raptor_search, query_wrong)
{
    cli_test_result result = execute_app("raptor", "search", "--query foo.fasta", "--index ", data("expected_results/b64_k19_w23_s8m.ibf"), "--output search.out");
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Validation failed for option --query: The file \"foo.fasta\" does not exist!\n"});
}

TEST_F(raptor_search, output_missing)
{
    cli_test_result result = execute_app("raptor", "search", "--query ", data("example_data/64/reads/all.fastq"), "--index ", data("expected_results/b64_k19_w23_s8m.ibf"));
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Option --output is required but not set.\n"});
}

TEST_F(raptor_search, output_wrong)
{
    cli_test_result result = execute_app("raptor", "search", "--query ", data("example_data/64/reads/all.fastq"), "--index ", data("expected_results/b64_k19_w23_s8m.ibf"), "--output foo/search.out");
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Validation failed for option --output: Cannot write \"foo/search.out\"!\n"});
}

TEST_F(raptor_search, kmer_window)
{
    cli_test_result result = execute_app("raptor", "search", "--kmer 20", "--window 19","--query ", data("example_data/64/reads/all.fastq"), "--index ", data("expected_results/b64_k19_w23_s8m.ibf"), "--output search.out");
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] The k-mer size cannot be bigger than the window size.\n"});
}

#include <seqan3/std/ranges>     // range comparisons
#include <string>                // strings
#include <vector>                // vectors

#include <seqan3/test/snippet/create_temporary_snippet_file.hpp>

seqan3::test::create_temporary_snippet_file tmp_ibf_file{"tmp.ibf", "dummy"};

#include "cli_test.hpp"

TEST_F(raptor, no_options)
{
    cli_test_result const result = execute_app("raptor");
    std::string const expected
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
    cli_test_result const result = execute_app("raptor", "build");
    std::string const expected
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
    cli_test_result const result = execute_app("raptor", "search");
    std::string const expected
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
    cli_test_result const result = execute_app("raptor", "foo");
    std::string const expected
    {
        "[Error] You either forgot or misspelled the subcommand! Please specify which sub-program you want to use: one "
        "of [build,search]. Use -h/--help for more information.\n"
    };
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(raptor, unknown_option)
{
    cli_test_result const result = execute_app("raptor", "-v");
    std::string const expected
    {
        "[Error] You either forgot or misspelled the subcommand! Please specify which sub-program you want to use: one "
        "of [build,search]. Use -h/--help for more information.\n"
    };
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(raptor_build, output_missing)
{
    cli_test_result const result = execute_app("raptor", "build",
                                                         "--size 8m",
                                                         data("bin1.fa"));
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Option --output is required but not set.\n"});
}

TEST_F(raptor_build, output_wrong)
{
    cli_test_result const result = execute_app("raptor", "build",
                                                         "--size 8m",
                                                         "--output foo/out.ibf",
                                                         data("bin1.fa"));
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Cannot write \"foo/out.ibf\"!\n"});
}

TEST_F(raptor_build, directory_missing)
{
    cli_test_result const result = execute_app("raptor", "build",
                                                         "--size 8m",
                                                         "--compute-minimiser",
                                                         data("bin1.fa"));
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Option --output is required but not set.\n"});
}

TEST_F(raptor_build, directory_wrong)
{
    cli_test_result const result = execute_app("raptor", "build",
                                                         "--size 8m",
                                                         "--compute-minimiser",
                                                         "--output foo/bar",
                                                         data("bin1.fa"));
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Cannot create directory: \"foo/bar\"!\n"});
}

TEST_F(raptor_build, size_missing)
{
    cli_test_result const result = execute_app("raptor", "build",
                                                         "--output ./ibf.out",
                                                         data("bin1.fa"));
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Option --size is required but not set.\n"});
}

TEST_F(raptor_build, size_wrong_space)
{
    cli_test_result const result = execute_app("raptor", "build",
                                                         "--size 8 m",
                                                         "--output ./ibf.out",
                                                         data("bin1.fa"));
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Validation failed for option --size: Value 8 must be an integer "
                                      "followed by [k,m,g,t] (case insensitive).\n"});
}

TEST_F(raptor_build, size_wrong_suffix)
{
    cli_test_result const result = execute_app("raptor", "build",
                                                         "--size 8x",
                                                         "--output ibf.out",
                                                         data("bin1.fa"));
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Validation failed for option --size: Value 8x must be an integer "
                                      "followed by [k,m,g,t] (case insensitive).\n"});
}

TEST_F(raptor_build, kmer_window)
{
    cli_test_result const result = execute_app("raptor", "build",
                                                         "--kmer 20",
                                                         "--window 19",
                                                         "--size 8m",
                                                         "--output ibf.out",
                                                         data("bin1.fa"));
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] The k-mer size cannot be bigger than the window size.\n"});
}

TEST_F(raptor_search, ibf_missing)
{
    cli_test_result const result = execute_app("raptor", "search",
                                                         "--query ", data("query.fq"),
                                                         "--output search.out");
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Option --index is required but not set.\n"});
}

TEST_F(raptor_search, ibf_wrong)
{
    cli_test_result const result = execute_app("raptor", "search",
                                                         "--query ", data("query.fq"),
                                                         "--index foo.ibf",
                                                         "--output search.out");
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Validation failed for option --index: The file \"foo.ibf\" does not "
                                      "exist!\n"});
}

TEST_F(raptor_search, query_missing)
{
    cli_test_result const result = execute_app("raptor", "search",
                                                         "--index ", tmp_ibf_file.file_path,
                                                         "--output search.out");
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Option --query is required but not set.\n"});
}

TEST_F(raptor_search, query_wrong)
{
    cli_test_result const result = execute_app("raptor", "search",
                                                         "--query foo.fasta",
                                                         "--index ", tmp_ibf_file.file_path,
                                                         "--output search.out");
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Validation failed for option --query: The file \"foo.fasta\" does not "
                                      "exist!\n"});
}

TEST_F(raptor_search, output_missing)
{
    cli_test_result const result = execute_app("raptor", "search",
                                                         "--query ", data("query.fq"),
                                                         "--index ", tmp_ibf_file.file_path);
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Option --output is required but not set.\n"});
}

TEST_F(raptor_search, output_wrong)
{
    cli_test_result const result = execute_app("raptor", "search",
                                                         "--query ", data("query.fq"),
                                                         "--index ", tmp_ibf_file.file_path,
                                                         "--output foo/search.out");
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Validation failed for option --output: Cannot write "
                                      "\"foo/search.out\"!\n"});
}

TEST_F(raptor_search, kmer_window)
{
    cli_test_result const result = execute_app("raptor", "search",
                                                         "--kmer 20",
                                                         "--window 19",
                                                         "--query ", data("query.fq"),
                                                         "--index ", tmp_ibf_file.file_path,
                                                         "--output search.out");
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] The k-mer size cannot be bigger than the window size.\n"});
}

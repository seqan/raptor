// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "inspect_test.hpp"

struct inspect_options : public inspect_test
{};

TEST_F(inspect_options, no_options)
{
    cli_test_result const result = execute_app("raptor", "inspect");
    std::string const expected{"Raptor-inspect - Prints the metadata of a Raptor index\n"
                               "======================================================\n"
                               "    raptor inspect --index <file> [--output <file>] [--format text|json]\n"
                               "    [--user-bins] [--detailed]\n"
                               "    Try -h or --help for more information.\n"};
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);
}

TEST_F(inspect_options, missing_index)
{
    cli_test_result const result = execute_app("raptor", "inspect", "--format json");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Option --index is required but not set.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(inspect_options, index_does_not_exist)
{
    cli_test_result const result = execute_app("raptor", "inspect", "--index does_not_exist.index");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] The file \"does_not_exist.index\" does not exist!\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(inspect_options, invalid_format)
{
    cli_test_result const result = execute_app("raptor", "inspect", "--index", ibf_path(16u, 19u), "--format xml");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err,
              std::string{"[Error] Validation failed for option --format: Value xml is not one of [text, json].\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(inspect_options, old_index)
{
    cli_test_result const result = execute_app("raptor", "inspect", "--index", data("3.0.index"));
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Unsupported index version. Check raptor upgrade.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(inspect_options, ambiguous_index)
{
    std::filesystem::copy_file(ibf_path(16u, 19u), "raptor.index");
    std::filesystem::copy_file(ibf_path(16u, 19u), "raptor.index_0");

    cli_test_result const result = execute_app("raptor", "inspect", "--index raptor.index");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err,
              std::string{"[Error] Ambiguous index. Both monolithic (raptor.index) and partitioned index "
                          "(raptor.index_0) exist. Please rename the monolithic index.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(inspect_options, output_file)
{
    cli_test_result const to_stdout = execute_app("raptor", "inspect", "--index", ibf_path(16u, 19u, is_hibf::yes));
    EXPECT_EQ(to_stdout.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(to_stdout);

    // The directory is created if it does not exist.
    cli_test_result const to_file = execute_app("raptor",
                                                "inspect",
                                                "--index",
                                                ibf_path(16u, 19u, is_hibf::yes),
                                                "--output directory/metadata.txt");
    EXPECT_EQ(to_file.out, std::string{});
    EXPECT_EQ(to_file.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(to_file);

    EXPECT_FALSE(to_stdout.out.empty());
    EXPECT_EQ(string_from_file("directory/metadata.txt"), to_stdout.out);

    // An existing file is overwritten.
    cli_test_result const overwrite = execute_app("raptor",
                                                  "inspect",
                                                  "--index",
                                                  ibf_path(16u, 19u, is_hibf::yes),
                                                  "--format json",
                                                  "--output directory/metadata.txt");
    RAPTOR_ASSERT_ZERO_EXIT(overwrite);
    EXPECT_TRUE(string_from_file("directory/metadata.txt").starts_with("{\n"));
}

// A failure leaves no output file behind.
TEST_F(inspect_options, no_output_on_error)
{
    std::ofstream{"not_an_index.index"} << "some content";

    cli_test_result const result =
        execute_app("raptor", "inspect", "--index not_an_index.index", "--output metadata.txt");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Unsupported index version. Check raptor upgrade.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
    EXPECT_FALSE(std::filesystem::exists("metadata.txt"));
}

// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "update_test.hpp"

struct update_options : public update_test
{};

// Without a subcommand, sharg prints a short description and exits successfully.
TEST_F(update_options, no_subcommand)
{
    cli_test_result const result = execute_app("raptor", "update");
    EXPECT_TRUE(result.out.starts_with("Raptor-update - Updates a Raptor index")) << result.out;
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);
}

TEST_F(update_options, unknown_subcommand)
{
    cli_test_result const result = execute_app("raptor", "update", "modify");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_TRUE(result.err.starts_with("[Error]")) << result.err;
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(update_options, insert_missing_index)
{
    cli_test_result const result =
        execute_app("raptor", "update", "insert", "--output out.index", "--insert", data("bin1.fa"));
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Option --index is required but not set.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(update_options, insert_missing_output)
{
    ASSERT_NO_FATAL_FAILURE(build_index(write_bin_file("bins.txt", bins(2u)), "raptor.index"));

    cli_test_result const result =
        execute_app("raptor", "update", "insert", "--index raptor.index", "--insert", data("bin1.fa"));
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Option --output is required but not set.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(update_options, insert_missing_insert)
{
    ASSERT_NO_FATAL_FAILURE(build_index(write_bin_file("bins.txt", bins(2u)), "raptor.index"));

    cli_test_result const result =
        execute_app("raptor", "update", "insert", "--index raptor.index", "--output out.index");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Option --insert is required but not set.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(update_options, delete_missing_delete)
{
    ASSERT_NO_FATAL_FAILURE(build_index(write_bin_file("bins.txt", bins(2u)), "raptor.index"));

    cli_test_result const result =
        execute_app("raptor", "update", "delete", "--index raptor.index", "--output out.index");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Option --delete is required but not set.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(update_options, index_does_not_exist)
{
    cli_test_result const result = execute_app("raptor",
                                               "update",
                                               "insert",
                                               "--index does_not_exist.index",
                                               "--output out.index",
                                               "--insert",
                                               data("bin1.fa"));
    EXPECT_EQ(result.out, std::string{});
    EXPECT_TRUE(result.err.starts_with("[Error] Validation failed for option --index:")) << result.err;
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(update_options, output_already_exists)
{
    ASSERT_NO_FATAL_FAILURE(build_index(write_bin_file("bins.txt", bins(2u)), "raptor.index"));
    std::ofstream{"out.index"} << "not empty";

    cli_test_result const result = execute_app("raptor",
                                               "update",
                                               "insert",
                                               "--index raptor.index",
                                               "--output out.index",
                                               "--insert",
                                               data("bin1.fa"));
    EXPECT_EQ(result.out, std::string{});
    EXPECT_TRUE(result.err.starts_with("[Error] Validation failed for option --output:")) << result.err;
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(update_options, empty_insert_file)
{
    ASSERT_NO_FATAL_FAILURE(build_index(write_bin_file("bins.txt", bins(2u)), "raptor.index"));
    std::ofstream{"empty.txt"};

    cli_test_result const result =
        execute_app("raptor", "update", "insert", "--index raptor.index", "--output out.index", "--insert empty.txt");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] The input file is empty.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

// Updating relies on the occupancy of the technical bins, which seqan::hibf only maintains if the layout reserved
// empty bins. The shipped indices are built without them.
TEST_F(update_options, index_without_empty_bins)
{
    cli_test_result const insertion = execute_app("raptor",
                                                  "update",
                                                  "insert",
                                                  "--index",
                                                  ibf_path(0u, 19u, is_hibf::yes),
                                                  "--output out.index",
                                                  "--insert",
                                                  data("bin1.fa"));
    EXPECT_EQ(insertion.out, std::string{});
    EXPECT_TRUE(insertion.err.starts_with("[Error] The index does not track the occupancy")) << insertion.err;
    RAPTOR_ASSERT_FAIL_EXIT(insertion);

    cli_test_result const deletion = execute_app("raptor",
                                                 "update",
                                                 "delete",
                                                 "--index",
                                                 ibf_path(0u, 19u, is_hibf::yes),
                                                 "--output out.index",
                                                 "--delete 0");
    EXPECT_EQ(deletion.out, std::string{});
    EXPECT_TRUE(deletion.err.starts_with("[Error] The index does not track the occupancy")) << deletion.err;
    RAPTOR_ASSERT_FAIL_EXIT(deletion);
}

TEST_F(update_options, ibf_index)
{
    cli_test_result const result = execute_app("raptor",
                                               "update",
                                               "insert",
                                               "--index",
                                               ibf_path(0u, 19u),
                                               "--output out.index",
                                               "--insert",
                                               data("bin1.fa"));
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Only an HIBF can be updated.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(update_options, multiple_files_per_user_bin)
{
    ASSERT_NO_FATAL_FAILURE(build_index(write_bin_file("bins.txt", bins(2u)), "raptor.index"));
    std::ofstream{"insert.txt"} << data("bin3.fa").string() << ' ' << data("bin4.fa").string() << '\n';

    cli_test_result const result =
        execute_app("raptor", "update", "insert", "--index raptor.index", "--output out.index", "--insert insert.txt");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_TRUE(result.err.find("multiple files per UB") != std::string::npos) << result.err;
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

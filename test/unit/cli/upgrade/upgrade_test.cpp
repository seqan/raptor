// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <raptor/test/cli_test.hpp>

struct upgrade : public raptor_base
{};

TEST_F(upgrade, via_fpr)
{
    cli_test_result const result =
        execute_app("raptor", "upgrade", "--input ", data("2.0.index"), "--output raptor.index", "--fpr 0.05");
    EXPECT_NE(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);

    compare_index(ibf_path(16, 19), "raptor.index");
}

TEST_F(upgrade, via_bin_path)
{
    std::filesystem::copy_file(data("bin1.fa"), std::filesystem::current_path() / "bin1.fa");
    std::filesystem::copy_file(data("bin2.fa"), std::filesystem::current_path() / "bin2.fa");
    std::filesystem::copy_file(data("bin3.fa"), std::filesystem::current_path() / "bin3.fa");
    std::filesystem::copy_file(data("bin4.fa"), std::filesystem::current_path() / "bin4.fa");

    cli_test_result const result =
        execute_app("raptor", "upgrade", "--input ", data("2.0.index"), "--output raptor.index");
    EXPECT_NE(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);

    compare_index(ibf_path(16, 19), "raptor.index");
}

TEST_F(upgrade, via_bin_file)
{
    { // generate input file
        std::ofstream file{"raptor_cli_test.txt"};
        for (auto && file_path : get_repeated_bins(16u))
            file << file_path << '\n';
        file << '\n';
    }

    cli_test_result const result = execute_app("raptor",
                                               "upgrade",
                                               "--input ",
                                               data("2.0.index"),
                                               "--output raptor.index",
                                               "--bins raptor_cli_test.txt");
    EXPECT_NE(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);

    compare_index(ibf_path(16, 19), "raptor.index");
}

TEST_F(upgrade, compressed)
{
    cli_test_result const result = execute_app("raptor",
                                               "upgrade",
                                               "--input ",
                                               data("2.0.compressed.index"),
                                               "--output raptor.index",
                                               "--fpr 0.05");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Compressed upgrade not yet supported on main branch.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(upgrade, partitioned_ibf_via_fpr)
{
    cli_test_result const result1 = execute_app("raptor",
                                                "upgrade",
                                                "--input ",
                                                data("2.0.partitioned.index"),
                                                "--output raptor.index",
                                                "--fpr 0.05");

    EXPECT_NE(result1.out, std::string{});
    EXPECT_EQ(result1.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result1);

    // raptor 3.0 has variable size partitions, so we cannot compare the indices
    cli_test_result const result2 = execute_app("raptor",
                                                "search",
                                                "--output search.out",
                                                "--error 1",
                                                "--index raptor.index",
                                                "--quiet",
                                                "--query ",
                                                data("query.fq"));
    EXPECT_EQ(result2.out, std::string{});
    EXPECT_EQ(result2.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result2);

    compare_search(16, 1, "search.out");
}

TEST_F(upgrade, partitioned_ibf_via_bin_path)
{
    std::filesystem::copy_file(data("bin1.fa"), std::filesystem::current_path() / "bin1.fa");
    std::filesystem::copy_file(data("bin2.fa"), std::filesystem::current_path() / "bin2.fa");
    std::filesystem::copy_file(data("bin3.fa"), std::filesystem::current_path() / "bin3.fa");
    std::filesystem::copy_file(data("bin4.fa"), std::filesystem::current_path() / "bin4.fa");

    cli_test_result const result1 =
        execute_app("raptor", "upgrade", "--input ", data("2.0.partitioned.index"), "--output raptor.index");

    EXPECT_NE(result1.out, std::string{});
    EXPECT_EQ(result1.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result1);

    // raptor 3.0 has variable size partitions, so we cannot compare the indices
    cli_test_result const result2 = execute_app("raptor",
                                                "search",
                                                "--output search.out",
                                                "--error 1",
                                                "--index raptor.index",
                                                "--quiet",
                                                "--query ",
                                                data("query.fq"));
    EXPECT_EQ(result2.out, std::string{});
    EXPECT_EQ(result2.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result2);

    compare_search(16, 1, "search.out");
}

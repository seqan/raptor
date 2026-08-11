// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "update_test.hpp"

struct update_insert : public update_test
{};

TEST_F(update_insert, single_sequence_file)
{
    ASSERT_NO_FATAL_FAILURE(build_index(write_bin_file("bins.txt", bins(3u)), "raptor.index"));

    cli_test_result const result = execute_app("raptor",
                                               "update",
                                               "insert",
                                               "--threads 1",
                                               "--index raptor.index",
                                               "--output updated.index",
                                               "--insert",
                                               data("bin4.fa"));
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);

    search_result const searched = search("updated.index");
    ASSERT_EQ(searched.user_bins.size(), 4u);
    EXPECT_TRUE(searched.user_bins.back().ends_with("bin4.fa"));
    expect_hits(searched, all_ids(4u));
}

TEST_F(update_insert, file_of_files)
{
    ASSERT_NO_FATAL_FAILURE(build_index(write_bin_file("bins.txt", bins(2u)), "raptor.index"));
    write_bin_file("insert.txt", {data("bin3.fa").string(), data("bin4.fa").string()});

    cli_test_result const result = execute_app("raptor",
                                               "update",
                                               "insert",
                                               "--threads 1",
                                               "--index raptor.index",
                                               "--output updated.index",
                                               "--insert insert.txt");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);

    expect_hits(search("updated.index"), all_ids(4u));
}

TEST_F(update_insert, input_index_is_not_modified)
{
    ASSERT_NO_FATAL_FAILURE(build_index(write_bin_file("bins.txt", bins(3u)), "raptor.index"));
    std::string const before = string_from_file("raptor.index", std::ios::binary);

    cli_test_result const result = execute_app("raptor",
                                               "update",
                                               "insert",
                                               "--threads 1",
                                               "--index raptor.index",
                                               "--output updated.index",
                                               "--insert",
                                               data("bin4.fa"));
    RAPTOR_ASSERT_ZERO_EXIT(result);

    EXPECT_EQ(before, string_from_file("raptor.index", std::ios::binary));
    expect_hits(search("raptor.index"), all_ids(3u));
}

TEST_F(update_insert, repeated_insertions)
{
    ASSERT_NO_FATAL_FAILURE(build_index(write_bin_file("bins.txt", bins(1u)), "raptor.index"));

    std::string current{"raptor.index"};
    for (size_t i = 2; i <= 4u; ++i)
    {
        std::string const next = "raptor." + std::to_string(i) + ".index";
        cli_test_result const result = execute_app("raptor",
                                                   "update",
                                                   "insert",
                                                   "--threads 1",
                                                   "--index",
                                                   current,
                                                   "--output",
                                                   next,
                                                   "--insert",
                                                   data("bin" + std::to_string(i) + ".fa"));
        EXPECT_EQ(result.err, std::string{}) << "Insertion " << i;
        RAPTOR_ASSERT_ZERO_EXIT(result);

        search_result const searched = search(next);
        ASSERT_EQ(searched.user_bins.size(), i) << "Insertion " << i;
        expect_hits(searched, all_ids(i));

        current = next;
    }
}

// Inserting as many user bins as the index already holds makes it outgrow its layout, which triggers rebuilds.
// With these parameters the index grows from one IBF to four.
TEST_F(update_insert, triggers_rebuild)
{
    std::vector<std::string> const initial = get_repeated_bins(4u) | std::ranges::to<std::vector<std::string>>();
    ASSERT_NO_FATAL_FAILURE(build_index(write_bin_file("bins.txt", initial), "raptor.index"));
    write_bin_file("insert.txt", initial); // The same files again, as new user bins.

    cli_test_result const result = execute_app("raptor",
                                               "update",
                                               "insert",
                                               "--threads 1",
                                               "--index raptor.index",
                                               "--output updated.index",
                                               "--insert insert.txt");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);

    search_result const searched = search("updated.index");
    ASSERT_EQ(searched.user_bins.size(), 2u * initial.size());
    expect_hits(searched, all_ids(2u * initial.size()));
}

TEST_F(update_insert, no_partial_rebuild)
{
    std::vector<std::string> const initial = get_repeated_bins(4u) | std::ranges::to<std::vector<std::string>>();
    ASSERT_NO_FATAL_FAILURE(build_index(write_bin_file("bins.txt", initial), "raptor.index"));
    write_bin_file("insert.txt", initial);

    cli_test_result const result = execute_app("raptor",
                                               "update",
                                               "insert",
                                               "--threads 1",
                                               "--index raptor.index",
                                               "--output updated.index",
                                               "--no-partial-rebuild",
                                               "--insert insert.txt");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);

    search_result const searched = search("updated.index");
    ASSERT_EQ(searched.user_bins.size(), 2u * initial.size());
    expect_hits(searched, all_ids(2u * initial.size()));
}

// The index and the files to insert may be of different types, hence the file type is deduced per user bin.
TEST_F(update_insert, minimiser_input)
{
    write_bin_file("bins.txt", bins(3u));
    write_bin_file("insert_source.txt", {data("bin4.fa").string()});

    for (std::string_view const name : {"bins.txt", "insert_source.txt"})
    {
        cli_test_result const prepare = execute_app("raptor",
                                                    "prepare",
                                                    "--threads 1",
                                                    "--quiet",
                                                    "--kmer 19",
                                                    "--window 19",
                                                    "--kmer-count-cutoff 1",
                                                    "--input",
                                                    name,
                                                    "--output",
                                                    std::string{"prepared_"} + std::string{name});
        ASSERT_EQ(prepare.err, std::string{});
        RAPTOR_ASSERT_ZERO_EXIT(prepare);
    }

    ASSERT_NO_FATAL_FAILURE(build_index("prepared_bins.txt/minimiser.list", "raptor.index"));

    cli_test_result const result = execute_app("raptor",
                                               "update",
                                               "insert",
                                               "--threads 1",
                                               "--index raptor.index",
                                               "--output updated.index",
                                               "--insert prepared_insert_source.txt/minimiser.list");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);

    search_result const searched = search("updated.index");
    ASSERT_EQ(searched.user_bins.size(), 4u);
    EXPECT_TRUE(searched.user_bins.back().ends_with("bin4.minimiser"));
    expect_hits(searched, all_ids(4u));
}

// Without errors, the queries do not match bin4.fa. The inserted bin must not turn up regardless.
TEST_F(update_insert, does_not_report_everything)
{
    ASSERT_NO_FATAL_FAILURE(build_index(write_bin_file("bins.txt", bins(3u)), "raptor.index"));

    cli_test_result const result = execute_app("raptor",
                                               "update",
                                               "insert",
                                               "--threads 1",
                                               "--index raptor.index",
                                               "--output updated.index",
                                               "--insert",
                                               data("bin4.fa"));
    RAPTOR_ASSERT_ZERO_EXIT(result);

    search_result const searched = search("updated.index", 0u);
    ASSERT_EQ(searched.user_bins.size(), 4u);
    expect_hits(searched, std::set<size_t>{0u, 1u, 2u});
}

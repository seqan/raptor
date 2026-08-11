// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "update_test.hpp"

struct update_delete : public update_test
{};

TEST_F(update_delete, single_user_bin)
{
    ASSERT_NO_FATAL_FAILURE(build_index(write_bin_file("bins.txt", bins(4u)), "raptor.index"));

    cli_test_result const result = execute_app("raptor",
                                               "update",
                                               "delete",
                                               "--threads 1",
                                               "--index raptor.index",
                                               "--output updated.index",
                                               "--delete 1");
    EXPECT_EQ(result.out, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);

    // The user bin IDs of the remaining bins do not change.
    expect_hits(search("updated.index"), std::set<size_t>{0u, 2u, 3u});
}

TEST_F(update_delete, multiple_user_bins)
{
    ASSERT_NO_FATAL_FAILURE(build_index(write_bin_file("bins.txt", bins(4u)), "raptor.index"));

    cli_test_result const result = execute_app("raptor",
                                               "update",
                                               "delete",
                                               "--threads 1",
                                               "--index raptor.index",
                                               "--output updated.index",
                                               "--delete 0",
                                               "--delete 2");
    EXPECT_EQ(result.out, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);

    expect_hits(search("updated.index"), std::set<size_t>{1u, 3u});
}

TEST_F(update_delete, input_index_is_not_modified)
{
    ASSERT_NO_FATAL_FAILURE(build_index(write_bin_file("bins.txt", bins(4u)), "raptor.index"));
    std::string const before = string_from_file("raptor.index", std::ios::binary);

    cli_test_result const result = execute_app("raptor",
                                               "update",
                                               "delete",
                                               "--threads 1",
                                               "--index raptor.index",
                                               "--output updated.index",
                                               "--delete 1");
    RAPTOR_ASSERT_ZERO_EXIT(result);

    EXPECT_EQ(before, string_from_file("raptor.index", std::ios::binary));
    expect_hits(search("raptor.index"), all_ids(4u));
}

/*!\brief The technical bins of a deleted user bin are reused by a later insertion.
 * \details
 * Both indices end up with the same six user bin paths, so their sizes are comparable. The one that could reuse the
 * freed technical bins must be the smaller one.
 */
TEST_F(update_delete, deleted_bins_are_reused)
{
    ASSERT_NO_FATAL_FAILURE(build_index(write_bin_file("bins.txt", bins(4u)), "raptor.index"));
    write_bin_file("insert.txt", {data("bin1.fa").string(), data("bin2.fa").string()});

    cli_test_result const deletion = execute_app("raptor",
                                                 "update",
                                                 "delete",
                                                 "--threads 1",
                                                 "--index raptor.index",
                                                 "--output deleted.index",
                                                 "--delete 0",
                                                 "--delete 1");
    RAPTOR_ASSERT_ZERO_EXIT(deletion);

    cli_test_result const with_reuse = execute_app("raptor",
                                                   "update",
                                                   "insert",
                                                   "--threads 1",
                                                   "--index deleted.index",
                                                   "--output with_reuse.index",
                                                   "--insert insert.txt");
    EXPECT_EQ(with_reuse.out, std::string{});
    EXPECT_EQ(with_reuse.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(with_reuse);

    cli_test_result const without_reuse = execute_app("raptor",
                                                      "update",
                                                      "insert",
                                                      "--threads 1",
                                                      "--index raptor.index",
                                                      "--output without_reuse.index",
                                                      "--insert insert.txt");
    RAPTOR_ASSERT_ZERO_EXIT(without_reuse);

    search_result const searched = search("with_reuse.index");
    ASSERT_EQ(searched.user_bins.size(), 6u);
    // The deleted user bins keep their IDs and stay empty; the inserted ones are appended.
    expect_hits(searched, std::set<size_t>{2u, 3u, 4u, 5u});

    EXPECT_LT(std::filesystem::file_size("with_reuse.index"), std::filesystem::file_size("without_reuse.index"))
        << "The freed technical bins were not reused.";
}

// Emptying an IBF completely retires it and every IBF below it. With these parameters the index holds four IBFs,
// three of which are retired. The index has to stay usable afterwards.
TEST_F(update_delete, all_user_bins)
{
    std::vector<std::string> const initial = get_repeated_bins(4u) | std::ranges::to<std::vector<std::string>>();
    ASSERT_NO_FATAL_FAILURE(build_index(write_bin_file("bins.txt", initial), "raptor.index"));
    write_bin_file("insert.txt", initial);

    // Insert to obtain an index with more than one IBF.
    cli_test_result const insertion = execute_app("raptor",
                                                  "update",
                                                  "insert",
                                                  "--threads 1",
                                                  "--index raptor.index",
                                                  "--output inserted.index",
                                                  "--insert insert.txt");
    ASSERT_EQ(insertion.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(insertion);

    size_t const number_of_user_bins = 2u * initial.size();
    std::string delete_arguments{};
    std::vector<size_t> deleted{};
    for (size_t i = 0; i < number_of_user_bins; ++i)
    {
        delete_arguments += "--delete " + std::to_string(i) + ' ';
        deleted.push_back(i);
    }

    cli_test_result const deletion = execute_app("raptor",
                                                 "update",
                                                 "delete",
                                                 "--threads 1",
                                                 "--index inserted.index",
                                                 "--output deleted.index",
                                                 delete_arguments);
    EXPECT_EQ(deletion.out, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(deletion);

    // Every query is still reported, but without any user bin.
    expect_hits(search("deleted.index"), std::set<size_t>{});

    // Inserting into an index whose IBFs have been retired must work.
    cli_test_result const reinsertion = execute_app("raptor",
                                                    "update",
                                                    "insert",
                                                    "--threads 1",
                                                    "--index deleted.index",
                                                    "--output updated.index",
                                                    "--insert",
                                                    data("bin1.fa"));
    EXPECT_EQ(reinsertion.out, std::string{});
    EXPECT_EQ(reinsertion.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(reinsertion);

    search_result const searched = search("updated.index");
    ASSERT_EQ(searched.user_bins.size(), number_of_user_bins + 1u);
    expect_hits(searched, std::set<size_t>{number_of_user_bins});
}

TEST_F(update_delete, then_insert)
{
    ASSERT_NO_FATAL_FAILURE(build_index(write_bin_file("bins.txt", bins(4u)), "raptor.index"));

    cli_test_result const deletion = execute_app("raptor",
                                                 "update",
                                                 "delete",
                                                 "--threads 1",
                                                 "--index raptor.index",
                                                 "--output deleted.index",
                                                 "--delete 3");
    RAPTOR_ASSERT_ZERO_EXIT(deletion);

    cli_test_result const insertion = execute_app("raptor",
                                                  "update",
                                                  "insert",
                                                  "--threads 1",
                                                  "--index deleted.index",
                                                  "--output updated.index",
                                                  "--insert",
                                                  data("bin4.fa"));
    EXPECT_EQ(insertion.out, std::string{});
    EXPECT_EQ(insertion.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(insertion);

    search_result const searched = search("updated.index");
    ASSERT_EQ(searched.user_bins.size(), 5u);
    expect_hits(searched, std::set<size_t>{0u, 1u, 2u, 4u});
}

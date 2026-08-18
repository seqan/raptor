// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/io/sequence_file/input.hpp>

#include <raptor/dna4_traits.hpp>

#include "update_test.hpp"

struct update_rebuild : public update_test
{
    // Writes one 250 bp read per user bin, taken from the beginning of that bin's own sequence.
    static std::string write_queries(std::string const & filename, std::vector<std::string> const & bins)
    {
        static constexpr size_t query_length{250u};

        seqan3::sequence_file_output fout{filename}; // FASTQ, deduced from the .fq extension.
        std::vector<seqan3::phred42> const quality(query_length, seqan3::assign_rank_to(40u, seqan3::phred42{}));

        for (size_t i = 0; i < bins.size(); ++i)
        {
            seqan3::sequence_file_input<raptor::dna4_traits, seqan3::fields<seqan3::field::seq>> fin{bins[i]};
            auto & sequence = (*fin.begin()).sequence();
            EXPECT_GE(sequence.size(), query_length) << bins[i];

            std::vector<seqan3::dna4> const read(sequence.begin(), sequence.begin() + query_length);
            fout.emplace_back(read, "query_" + std::to_string(i), quality);
        }

        return filename;
    }

    // Like update_test::search, which is bound to query.fq.
    search_result search_queries(std::string const & index, std::string const & queries)
    {
        std::string const output = index + ".search";
        cli_test_result const result = execute_app("raptor",
                                                   "search",
                                                   "--quiet",
                                                   "--error 0",
                                                   "--query_length 250",
                                                   "--index",
                                                   index,
                                                   "--query",
                                                   queries,
                                                   "--output",
                                                   output);
        EXPECT_EQ(result.err, std::string{});
        EXPECT_EQ(result.exit_code, 0) << result.command;

        return parse_search_result(output);
    }

    // Checks that query `i` reports user bin `i`, and that no deleted user bin is reported anywhere.
    // The query file may hold more queries than the index holds user bins, e.g. while inserting batch by batch.
    static void expect_own_bin(search_result const & result,
                               size_t const number_of_user_bins,
                               std::set<size_t> const & deleted = {})
    {
        ASSERT_GE(result.hits.size(), number_of_user_bins);

        for (size_t i = 0; i < number_of_user_bins; ++i)
        {
            auto const it = result.hits.find("query_" + std::to_string(i));
            ASSERT_NE(it, result.hits.end()) << "Query " << i << " was not reported.";

            if (!deleted.contains(i))
            {
                EXPECT_TRUE(it->second.contains(i)) << "User bin " << i << " was not found by its own query.";
            }

            for (size_t const hit : it->second)
                EXPECT_FALSE(deleted.contains(hit)) << "Deleted user bin " << hit << " was reported.";
        }
    }
};

// Covers the insertion path through the lower levels: a user bin is written into an IBF below the top level and
// into every merged bin above it, and IBFs are grown to make room. The empty bin fraction of 0.3 is what makes the
// layout wide enough for the lower levels to be used.
TEST_F(update_rebuild, insert_into_multi_level_index)
{
    std::vector<std::string> const all = generate_bins("bins", 200u);
    std::vector<std::string> const initial{all.begin(), all.begin() + 100};
    std::vector<std::string> const inserted{all.begin() + 100, all.end()};

    ASSERT_NO_FATAL_FAILURE(build_index(write_bin_file("bins.txt", initial), "raptor.index", 64u, "0.3"));
    write_bin_file("insert.txt", inserted);
    std::string const queries = write_queries("queries.fq", all);

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

    search_result const searched = search_queries("updated.index", queries);
    ASSERT_EQ(searched.user_bins.size(), all.size());
    ASSERT_NO_FATAL_FAILURE(expect_own_bin(searched, all.size()));
}

// Note that this workload needs no repair at all, so the flag currently changes nothing about what runs. It guards
// against the flag breaking insertion into a multi-level index.
TEST_F(update_rebuild, insert_into_multi_level_index_no_partial_rebuild)
{
    std::vector<std::string> const all = generate_bins("bins", 200u);
    std::vector<std::string> const initial{all.begin(), all.begin() + 100};
    std::vector<std::string> const inserted{all.begin() + 100, all.end()};

    ASSERT_NO_FATAL_FAILURE(build_index(write_bin_file("bins.txt", initial), "raptor.index", 64u, "0.3"));
    write_bin_file("insert.txt", inserted);
    std::string const queries = write_queries("queries.fq", all);

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

    search_result const searched = search_queries("updated.index", queries);
    ASSERT_EQ(searched.user_bins.size(), all.size());
    ASSERT_NO_FATAL_FAILURE(expect_own_bin(searched, all.size()));
}

TEST_F(update_rebuild, repeated_insertions_into_multi_level_index)
{
    std::vector<std::string> const all = generate_bins("bins", 200u);
    std::vector<std::string> const initial{all.begin(), all.begin() + 80};

    ASSERT_NO_FATAL_FAILURE(build_index(write_bin_file("bins.txt", initial), "raptor.index", 64u, "0.3"));
    std::string const queries = write_queries("queries.fq", all);

    std::string current{"raptor.index"};
    for (size_t batch = 0; batch < 3u; ++batch)
    {
        size_t const first = 80u + batch * 40u;
        write_bin_file("insert.txt", {all.begin() + first, all.begin() + first + 40u});
        std::string const next = "raptor." + std::to_string(batch) + ".index";

        cli_test_result const result = execute_app("raptor",
                                                   "update",
                                                   "insert",
                                                   "--threads 1",
                                                   "--index",
                                                   current,
                                                   "--output",
                                                   next,
                                                   "--insert insert.txt");
        EXPECT_EQ(result.err, std::string{}) << "Batch " << batch;
        RAPTOR_ASSERT_ZERO_EXIT(result);

        search_result const searched = search_queries(next, queries);
        ASSERT_EQ(searched.user_bins.size(), first + 40u) << "Batch " << batch;
        ASSERT_NO_FATAL_FAILURE(expect_own_bin(searched, first + 40u)) << "Batch " << batch;

        current = next;
    }
}

// An IBF whose bins all became empty is retired together with the subtree below it, see detail::retire_ibf. The
// index has to stay usable afterwards.
TEST_F(update_rebuild, delete_all_user_bins_of_multi_level_index)
{
    std::vector<std::string> const all = generate_bins("bins", 120u);

    ASSERT_NO_FATAL_FAILURE(build_index(write_bin_file("bins.txt", all), "raptor.index", 64u, "0.3"));
    std::string const queries = write_queries("queries.fq", all);

    std::string delete_arguments{};
    for (size_t i = 0; i < all.size(); ++i)
        delete_arguments += "--delete " + std::to_string(i) + ' ';

    cli_test_result const deletion = execute_app("raptor",
                                                 "update",
                                                 "delete",
                                                 "--threads 1",
                                                 "--index raptor.index",
                                                 "--output deleted.index",
                                                 delete_arguments);
    EXPECT_EQ(deletion.out, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(deletion);

    // Every query is still reported, but no user bin is.
    for (auto const & [query, hits] : search_queries("deleted.index", queries).hits)
        EXPECT_TRUE(hits.empty()) << "Query: " << query;

    cli_test_result const insertion = execute_app("raptor",
                                                  "update",
                                                  "insert",
                                                  "--threads 1",
                                                  "--index deleted.index",
                                                  "--output updated.index",
                                                  "--insert",
                                                  all.front());
    EXPECT_EQ(insertion.out, std::string{});
    EXPECT_EQ(insertion.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(insertion);

    search_result const searched = search_queries("updated.index", queries);
    ASSERT_EQ(searched.user_bins.size(), all.size() + 1u);
    // The re-inserted bin holds the same sequence as user bin 0, so query 0 finds it under its new ID.
    EXPECT_TRUE(searched.hits.at("query_0").contains(all.size()));
}

TEST_F(update_rebuild, delete_and_insert_in_multi_level_index)
{
    std::vector<std::string> const all = generate_bins("bins", 120u);

    ASSERT_NO_FATAL_FAILURE(build_index(write_bin_file("bins.txt", all), "raptor.index", 64u, "0.3"));
    std::string const queries = write_queries("queries.fq", all);

    std::set<size_t> const deleted{10u, 11u, 12u, 40u, 41u, 90u};
    std::string delete_arguments{};
    for (size_t const user_bin : deleted)
        delete_arguments += "--delete " + std::to_string(user_bin) + ' ';

    cli_test_result const deletion = execute_app("raptor",
                                                 "update",
                                                 "delete",
                                                 "--threads 1",
                                                 "--index raptor.index",
                                                 "--output deleted.index",
                                                 delete_arguments);
    EXPECT_EQ(deletion.out, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(deletion);

    ASSERT_NO_FATAL_FAILURE(expect_own_bin(search_queries("deleted.index", queries), all.size(), deleted));

    write_bin_file("insert.txt", {all.begin(), all.begin() + 6});
    cli_test_result const insertion = execute_app("raptor",
                                                  "update",
                                                  "insert",
                                                  "--threads 1",
                                                  "--index deleted.index",
                                                  "--output updated.index",
                                                  "--insert insert.txt");
    EXPECT_EQ(insertion.out, std::string{});
    EXPECT_EQ(insertion.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(insertion);

    search_result const searched = search_queries("updated.index", queries);
    ASSERT_EQ(searched.user_bins.size(), all.size() + 6u);
    // The deleted user bins keep their IDs and stay empty.
    ASSERT_NO_FATAL_FAILURE(expect_own_bin(searched, all.size(), deleted));
}

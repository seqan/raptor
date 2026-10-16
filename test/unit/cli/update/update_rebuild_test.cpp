// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <fstream>

#include <cereal/archives/binary.hpp>

#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/io/sequence_file/input.hpp>

#include <raptor/dna4_traits.hpp>
#include <raptor/index.hpp>

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

    /*!\brief Checks that the index holds exactly the user bins of `bins` that are not in `deleted`, in their order.
     * \details
     * After a full rebuild, the index is the same as one built from the remaining user bins only. Query `i` has to
     * report the new ID of user bin `bins[i]`, and the query of a deleted user bin nothing at all.
     */
    static void expect_remaining(search_result const & result,
                                 std::vector<std::string> const & bins,
                                 std::set<size_t> const & deleted)
    {
        std::vector<std::string> remaining{};
        for (size_t i = 0; i < bins.size(); ++i)
            if (!deleted.contains(i))
                remaining.push_back(bins[i]);
        ASSERT_EQ(result.user_bins, remaining);

        size_t user_bin_id{};
        for (size_t i = 0; i < bins.size(); ++i)
        {
            auto const it = result.hits.find("query_" + std::to_string(i));
            ASSERT_NE(it, result.hits.end()) << "Query " << i << " was not reported.";

            if (deleted.contains(i))
            {
                EXPECT_TRUE(it->second.empty()) << "The query of deleted user bin " << i << " has hits.";
            }
            else
            {
                EXPECT_TRUE(it->second.contains(user_bin_id)) << "User bin " << i << " was not found by its own "
                                                              << "query under ID " << user_bin_id << '.';
                ++user_bin_id;
            }
        }
    }

    //!\brief Deletes `user_bins` from `index`, writing the result to `output`.
    void delete_user_bins(std::string const & index, std::string const & output, std::set<size_t> const & user_bins)
    {
        std::string delete_arguments{};
        for (size_t const user_bin : user_bins)
            delete_arguments += "--delete " + std::to_string(user_bin) + ' ';

        cli_test_result const result = execute_app("raptor",
                                                   "update",
                                                   "delete",
                                                   "--threads 1",
                                                   "--index",
                                                   index,
                                                   "--output",
                                                   output,
                                                   delete_arguments);
        EXPECT_EQ(result.out, std::string{});
        EXPECT_EQ(result.err, std::string{});
        RAPTOR_ASSERT_ZERO_EXIT(result);
    }

    //!\brief Loads the HIBF index stored in `filename`.
    static raptor::raptor_index<raptor::index_structure::hibf> load_hibf(std::string const & filename)
    {
        raptor::raptor_index<raptor::index_structure::hibf> index{};
        std::ifstream stream{filename, std::ios::binary};
        cereal::BinaryInputArchive archive{stream};
        archive(index);
        return index;
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

// A full rebuild treats the index as if the deleted user bins never existed: it drops them from the bin path and
// renumbers the remaining user bins in their order. `--no-partial-rebuild` turns every repair into a full rebuild.
TEST_F(update_rebuild, full_rebuild_drops_deleted_user_bins)
{
    std::vector<std::string> const all = generate_bins("bins", 241u);
    std::vector<std::string> const initial{all.begin(), all.begin() + 120};

    ASSERT_NO_FATAL_FAILURE(build_index(write_bin_file("bins.txt", initial), "raptor.index"));
    std::string const queries = write_queries("queries.fq", all);

    std::set<size_t> const deleted{10u, 11u, 12u, 40u, 41u, 90u};
    ASSERT_NO_FATAL_FAILURE(delete_user_bins("raptor.index", "deleted.index", deleted));

    write_bin_file("insert.txt", {all.begin() + 120, all.begin() + 240});
    cli_test_result const insertion = execute_app("raptor",
                                                  "update",
                                                  "insert",
                                                  "--threads 1",
                                                  "--index deleted.index",
                                                  "--output updated.index",
                                                  "--no-partial-rebuild",
                                                  "--insert insert.txt");
    EXPECT_EQ(insertion.out, std::string{});
    EXPECT_EQ(insertion.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(insertion);

    // Only a full rebuild computes a new layout.
    size_t const layout_user_bins = load_hibf("updated.index").config().number_of_user_bins;
    ASSERT_NE(layout_user_bins, initial.size()) << "No full rebuild happened.";
    EXPECT_EQ(layout_user_bins, 240u - deleted.size());

    ASSERT_NO_FATAL_FAILURE(
        expect_remaining(search_queries("updated.index", queries), {all.begin(), all.begin() + 240}, deleted));

    // The next insertion appends to the renumbered user bins.
    cli_test_result const reinsertion = execute_app("raptor",
                                                    "update",
                                                    "insert",
                                                    "--threads 1",
                                                    "--index updated.index",
                                                    "--output reinserted.index",
                                                    "--insert",
                                                    all.back());
    EXPECT_EQ(reinsertion.out, std::string{});
    EXPECT_EQ(reinsertion.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(reinsertion);

    ASSERT_NO_FATAL_FAILURE(expect_remaining(search_queries("reinserted.index", queries), all, deleted));
}

// Like full_rebuild_drops_deleted_user_bins, but the full rebuild follows the insertion of a user bin: a single user
// bin far larger than the others still finds a place in the index, but the index has to be rebuilt in full
// afterwards. The inserted user bin must be part of the rebuild, the deleted ones must not.
TEST_F(update_rebuild, full_rebuild_after_insertion_drops_deleted_user_bins)
{
    std::vector<std::string> all = generate_bins("bins", 120u);
    std::string const large = generate_bins("large", 1u, 7u, {100'000u}).front();
    all.push_back(large);

    ASSERT_NO_FATAL_FAILURE(
        build_index(write_bin_file("bins.txt", {all.begin(), all.end() - 1}), "raptor.index", 64u, "0.3"));
    std::string const queries = write_queries("queries.fq", all);

    std::set<size_t> const deleted{3u, 50u, 77u};
    ASSERT_NO_FATAL_FAILURE(delete_user_bins("raptor.index", "deleted.index", deleted));

    cli_test_result const insertion = execute_app("raptor",
                                                  "update",
                                                  "insert",
                                                  "--threads 1",
                                                  "--index deleted.index",
                                                  "--output updated.index",
                                                  "--insert",
                                                  large);
    EXPECT_EQ(insertion.out, std::string{});
    EXPECT_EQ(insertion.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(insertion);

    // Only a full rebuild computes a new layout.
    size_t const layout_user_bins = load_hibf("updated.index").config().number_of_user_bins;
    ASSERT_NE(layout_user_bins, all.size() - 1u) << "No full rebuild happened.";
    EXPECT_EQ(layout_user_bins, all.size() - deleted.size());

    ASSERT_NO_FATAL_FAILURE(expect_remaining(search_queries("updated.index", queries), all, deleted));
}

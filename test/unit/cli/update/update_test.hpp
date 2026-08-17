// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides a fixture for `raptor update` CLI tests.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <array>
#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sequence_file/output.hpp>

#include <raptor/test/cli_test.hpp>

/*!\brief Fixture for `raptor update` tests.
 * \details
 * `raptor update` requires an index whose layout reserved empty bins, so the shipped indices cannot be used and each
 * test builds its own. The bins are cheap: `bin1.fa` to `bin4.fa` are a few kibibytes each.
 *
 * Searching `query.fq` with two errors reports every user bin of the index, and with zero errors every user bin
 * except the ones holding `bin4.fa`. Both are used as an oracle: after an update, the reported user bins must be
 * exactly the ones the index is expected to contain.
 */
struct update_test : public raptor_base
{
    //!\brief The result of searching an index, as user bin IDs rather than paths.
    struct search_result
    {
        //!\brief Maps a user bin ID to the file it holds.
        std::vector<std::string> user_bins{};
        //!\brief The user bins reported for each query, keyed by query name.
        std::map<std::string, std::set<size_t>> hits{};
    };

    //!\brief Writes `paths` to `filename`, one per line, and returns `filename`.
    static std::string write_bin_file(std::string const & filename, std::vector<std::string> const & paths)
    {
        std::ofstream file{filename};
        for (std::string const & path : paths)
            file << path << '\n';
        return filename;
    }

    //!\brief The paths of `bin1.fa` to `bin<count>.fa`.
    static std::vector<std::string> bins(size_t const count)
    {
        std::vector<std::string> result{};
        for (size_t i = 1; i <= count; ++i)
            result.push_back(cli_test::data("bin" + std::to_string(i) + ".fa").string());
        return result;
    }

    /*!\brief Writes `count` FASTA files with deterministic content into `directory` and returns their paths.
     * \details
     * The four shipped bins are identical whenever repeated, so a merged bin covering many of them holds barely more
     * distinct k-mers than a single one. Growing an index to several IBFs, and filling those IBFs, needs user bins
     * that actually differ. `std::mt19937_64` is fully specified, so the files are the same on every platform.
     */
    static std::vector<std::string>
    generate_bins(std::string const & directory, size_t const count, uint64_t const seed = 42u)
    {
        static constexpr std::array<size_t, 4u> lengths{500u, 1000u, 2000u, 4000u};

        std::filesystem::create_directories(directory);
        std::mt19937_64 rng{seed};
        std::vector<std::string> paths{};
        paths.reserve(count);

        for (size_t i = 0; i < count; ++i)
        {
            std::string const id{"bin_" + std::to_string(i)};
            std::string const path{directory + '/' + id + ".fa"};

            std::vector<seqan3::dna4> sequence(lengths[rng() % lengths.size()]);
            for (seqan3::dna4 & base : sequence)
                base.assign_rank(rng() % seqan3::dna4::alphabet_size);

            seqan3::sequence_file_output fout{path};
            fout.emplace_back(sequence, id);

            paths.push_back(path);
        }

        return paths;
    }

    void build_index(std::string const & bin_file,
                     std::string const & output,
                     size_t const tmax = 64u,
                     std::string const & empty_bin_fraction = "0.1")
    {
        cli_test_result const layout = execute_app("raptor",
                                                   "layout",
                                                   "--threads 1",
                                                   "--kmer 19",
                                                   "--window 19",
                                                   "--fpr 0.05",
                                                   "--hash 2",
                                                   "--disable-estimate-union",
                                                   "--empty-bin-fraction",
                                                   empty_bin_fraction,
                                                   "--tmax",
                                                   std::to_string(tmax),
                                                   "--input",
                                                   bin_file,
                                                   "--output",
                                                   output + ".layout");
        ASSERT_EQ(layout.err, std::string{});
        RAPTOR_ASSERT_ZERO_EXIT(layout);

        cli_test_result const build =
            execute_app("raptor", "build", "--threads 1", "--quiet", "--input", output + ".layout", "--output", output);
        ASSERT_EQ(build.err, std::string{});
        RAPTOR_ASSERT_ZERO_EXIT(build);
    }

    //!\brief Searches `query.fq` in `index` and parses the result.
    search_result search(std::string const & index, size_t const errors = 2u)
    {
        std::string const output = index + ".search";
        cli_test_result const result = execute_app("raptor",
                                                   "search",
                                                   "--quiet",
                                                   "--error",
                                                   std::to_string(errors),
                                                   "--index",
                                                   index,
                                                   "--query",
                                                   cli_test::data("query.fq"),
                                                   "--output",
                                                   output);
        EXPECT_EQ(result.err, std::string{});
        EXPECT_EQ(result.exit_code, 0) << result.command;

        return parse_search_result(output);
    }

    //!\brief Parses the output of `raptor search`.
    static search_result parse_search_result(std::string const & filename)
    {
        search_result result{};
        std::ifstream file{filename};
        std::string line{};

        while (std::getline(file, line))
        {
            if (line.starts_with("##")) // Parameter information.
                continue;

            if (line == "#QUERY_NAME\tUSER_BINS")
                break;

            // A user bin: `#<id>\t<path>`
            size_t const tab = line.find('\t');
            EXPECT_NE(tab, std::string::npos) << line;
            EXPECT_EQ(std::stoull(line.substr(1u, tab - 1u)), result.user_bins.size()) << "User bin IDs are not dense.";
            result.user_bins.push_back(line.substr(tab + 1u));
        }

        while (std::getline(file, line))
        {
            size_t const tab = line.find('\t');
            EXPECT_NE(tab, std::string::npos) << line;
            std::string const query = line.substr(0u, tab);
            std::set<size_t> & hits = result.hits[query];

            std::string const ids = line.substr(tab + 1u);
            for (size_t start = 0u; start < ids.size();)
            {
                size_t const comma = std::min(ids.find(',', start), ids.size());
                hits.emplace(std::stoull(ids.substr(start, comma - start)));
                start = comma + 1u;
            }
        }

        return result;
    }

    //!\brief The user bin IDs `[0, count)`.
    static std::set<size_t> all_ids(size_t const count)
    {
        std::set<size_t> result{};
        for (size_t i = 0; i < count; ++i)
            result.emplace(i);
        return result;
    }

    //!\brief Checks that every query reports exactly `expected` and that the index holds `expected.size()` user bins.
    static void expect_hits(search_result const & result, std::set<size_t> const & expected)
    {
        EXPECT_FALSE(result.hits.empty()) << "No queries were reported.";
        for (auto const & [query, hits] : result.hits)
            EXPECT_EQ(hits, expected) << "Query: " << query;
    }
};

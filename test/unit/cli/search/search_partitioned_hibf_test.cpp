// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include <seqan3/test/tmp_directory.hpp>

#include <raptor/test/cli_test.hpp>

struct search_partitioned_hibf :
    public raptor_base,
    public testing::WithParamInterface<std::tuple<size_t, size_t, std::string>>
{
    /*!\brief The length of each generated user bin.
     * \details
     * A partitioned layout computes MinHash sketches, which need 16 buckets of 40 values each. Bucketing on the low
     * bits of the seeded canonical k-mer value starves one bucket to about 6 % of its share, so a user bin needs
     * roughly 10'000 distinct 19-mers. 20'000 bases leave a margin of about 2x.
     */
    static constexpr size_t bin_length{20'000u};
    //!\brief The length of each query.
    static constexpr size_t query_length{100u};
    //!\brief The seed for the user bins; the queries without hits use another one.
    static constexpr uint64_t bin_seed{42u};

    //!\brief Used for all temporary files within a test.
    seqan3::test::tmp_directory const test_directory{};

    //!\brief The result of `raptor search`.
    struct search_result
    {
        //!\brief The path of each user bin, indexed by user bin ID.
        std::vector<std::string> user_bins{};
        //!\brief The name of each query, in the order of the search output.
        std::vector<std::string> queries{};
        //!\brief The sorted IDs of the user bins reported for each query; `hits[i]` belongs to `queries[i]`.
        std::vector<std::vector<size_t>> hits{};
    };

    //!\brief A random DNA sequence of `length` bases. `std::mt19937_64` is fully specified, so it is portable.
    static std::string random_sequence(std::mt19937_64 & rng, size_t const length)
    {
        std::string sequence(length, 'A');
        for (char & base : sequence)
            base = "ACGT"[rng() % 4u];
        return sequence;
    }

    //!\brief The user bins, generated once for all tests.
    static inline std::vector<std::string> const user_bins = []()
    {
        std::mt19937_64 rng{bin_seed};
        std::vector<std::string> result{};
        size_t const number_of_user_bins{128};

        for (size_t i = 0; i < number_of_user_bins; ++i)
        {
            result.push_back(random_sequence(rng, bin_length));
        }

        return result;
    }();

    //!\brief The queries for the user bins.
    struct query_sets
    {
        //!\brief `zero_errors[i]` is cut from user bin `i`.
        std::vector<std::string> zero_errors{};
        //!\brief `one_error[i]` is `zero_errors[i]` with one substitution.
        std::vector<std::string> one_error{};
        //!\brief Random sequences; none is expected to be found.
        std::vector<std::string> no_hits{};
    };

    /*!\brief The queries, generated once for all tests.
     * \details
     * Not a namespace-scope structured binding: referencing one from a function with OpenMP enabled crashes
     * Clang 19 and icpx 2025.0 in `SemaOpenMP::isOpenMPCapturedDecl`.
     */
    static inline query_sets const queries = []()
    {
        std::mt19937_64 rng{bin_seed + 1u};
        std::vector<std::string> result_0_errors{};
        std::vector<std::string> result_1_errors{};
        std::vector<std::string> result_no_hit{};

        for (size_t i = 0; i < user_bins.size(); ++i)
        {
            size_t const offset = rng() % (bin_length - query_length + 1u);
            std::string query = user_bins[i].substr(offset, query_length);

            result_0_errors.push_back(query);

            char & base = query[query_length / 2];
            base = (base == 'A') ? 'C' : 'A';

            result_1_errors.push_back(query);
            result_no_hit.push_back(random_sequence(rng, query_length));
        }

        return query_sets{.zero_errors = result_0_errors, .one_error = result_1_errors, .no_hits = result_no_hit};
    }();

    //!\brief Writes the user bin file for `raptor layout`, one path per line.
    void write_user_bins(std::filesystem::path const & filename, std::vector<std::string> const & bins)
    {
        std::ofstream layout_input_file{filename};
        for (size_t i{0}; i < bins.size(); ++i)
        {
            std::string const name = "user_bin_" + std::to_string(i);
            auto const path = test_directory.path() / (name + ".fa");
            std::ofstream file{path};
            file << '>' << name << '\n' << bins[i] << '\n';
            layout_input_file << path.string() << '\n'; // `operator<<` would quote the path.
        }
    }

    /*!\brief Writes one query per user bin, cut from that bin, with `errors` substitutions.
     * \details Query `query_<i>` is expected to be reported for user bin `i` only.
     */
    static void write_queries(std::filesystem::path const & filename, std::vector<std::string> const & queries)
    {
        std::ofstream file{filename};

        size_t i{0};
        for (auto const & query : queries)
            file << ">query_" << i++ << '\n' << query << '\n';
    }

    //!\brief Parses the output of `raptor search`.
    static search_result parse_search_result(std::filesystem::path const & filename)
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
            std::string const & query = result.queries.emplace_back(line.substr(0u, tab));
            std::vector<size_t> & hits = result.hits.emplace_back();

            std::string const ids = line.substr(tab + 1u);
            for (size_t start = 0u; start < ids.size();)
            {
                size_t const comma = std::min(ids.find(',', start), ids.size());
                size_t const id = std::stoull(ids.substr(start, comma - start));
                EXPECT_LT(id, result.user_bins.size()) << "Query: " << query;
                hits.push_back(id);
                start = comma + 1u;
            }
            std::ranges::sort(hits);
        }

        return result;
    }

    //!\brief Runs `raptor layout` and `raptor build` on `input_filename`.
    void build_partitioned_index(std::filesystem::path const & input_filename,
                                 std::filesystem::path const & index_filename,
                                 size_t const parts,
                                 std::string const & approach)
    {
        std::filesystem::path const layout_filename = test_directory.path() / "raptor_cli_test.layout";

        { // build layout
            cli_test_result const result = execute_app("raptor",
                                                       "layout",
                                                       "--partitioning-approach",
                                                       approach,
                                                       "--kmer 19",
                                                       "--threads 1",
                                                       "--input",
                                                       input_filename,
                                                       "--tmax 64",
                                                       "--fpr 0.05",
                                                       "--relaxed-fpr 0.5",
                                                       "--number-of-partitions",
                                                       std::to_string(parts),
                                                       "--output",
                                                       layout_filename);

            EXPECT_EQ(result.out, std::string{});
            EXPECT_EQ(result.err, std::string{});
            RAPTOR_ASSERT_ZERO_EXIT(result);
        }

        ASSERT_TRUE(std::filesystem::exists(layout_filename));

        { // build index
            cli_test_result const result =
                execute_app("raptor", "build", "--output", index_filename, "--quiet", "--input", layout_filename);

            EXPECT_EQ(result.out, std::string{});
            EXPECT_EQ(result.err, std::string{});
            RAPTOR_ASSERT_ZERO_EXIT(result);
        }
    }
};

TEST_P(search_partitioned_hibf, with_threshold)
{
    auto const [parts, number_of_errors, approach] = GetParam();

    std::filesystem::path const input_filename = test_directory.path() / "raptor_cli_test.txt";
    std::filesystem::path const query_filename = test_directory.path() / "query.fa";
    std::filesystem::path const index_filename = test_directory.path() / "raptor.index";
    std::filesystem::path const search_filename = test_directory.path() / "search.out";

    write_user_bins(input_filename, user_bins);

    if (number_of_errors == 0)
        write_queries(query_filename, queries.zero_errors);
    else
        write_queries(query_filename, queries.one_error);

    ASSERT_NO_FATAL_FAILURE(build_partitioned_index(input_filename, index_filename, parts, approach));

    { // check with search if index contains expected input
        cli_test_result const result = execute_app("raptor",
                                                   "search",
                                                   "--output",
                                                   search_filename,
                                                   "--error ",
                                                   std::to_string(number_of_errors),
                                                   "--index",
                                                   index_filename,
                                                   "--quiet",
                                                   "--query",
                                                   query_filename);
        EXPECT_EQ(result.out, std::string{});
        EXPECT_EQ(result.err, std::string{});
        RAPTOR_ASSERT_ZERO_EXIT(result);
    }

    ASSERT_TRUE(std::filesystem::exists(search_filename));
    search_result const result = parse_search_result(search_filename);

    EXPECT_EQ(result.user_bins.size(), user_bins.size());
    EXPECT_EQ(result.queries.size(), user_bins.size());

    for (size_t i = 0; i < result.queries.size(); ++i)
    {
        ASSERT_EQ(result.hits[i].size(), 1) << "There should be only 1 hit but there are/is " << result.hits[i].size();

        std::string const & user_bin = result.user_bins[result.hits[i][0]];
        std::string const & query = result.queries[i];

        std::string const start = test_directory.path().string() + "/user_bin_";
        ASSERT_TRUE(user_bin.starts_with(start)) << user_bin;
        ASSERT_TRUE(query.starts_with("query_")) << query;

        size_t const query_idx = std::stoull(query.substr(6u));
        size_t const user_bin_idx = std::stoull(user_bin.substr(start.size()));

        EXPECT_EQ(query_idx, user_bin_idx) << "query_<j> is cut from user bin j";
    }
}

TEST_P(search_partitioned_hibf, no_hits)
{
    auto const [parts, number_of_errors, approach] = GetParam();

    std::filesystem::path const input_filename = test_directory.path() / "raptor_cli_test.txt";
    std::filesystem::path const query_filename = test_directory.path() / "query_no_hits.fa";
    std::filesystem::path const index_filename = test_directory.path() / "raptor.index";
    std::filesystem::path const search_filename = test_directory.path() / "search.out";

    write_user_bins(input_filename, user_bins);
    write_queries(query_filename, queries.no_hits);

    ASSERT_NO_FATAL_FAILURE(build_partitioned_index(input_filename, index_filename, parts, approach));

    {
        cli_test_result const result = execute_app("raptor",
                                                   "search",
                                                   "--output",
                                                   search_filename,
                                                   "--error ",
                                                   std::to_string(number_of_errors),
                                                   "--tau 0.99",
                                                   "--index ",
                                                   index_filename,
                                                   "--quiet",
                                                   "--query ",
                                                   query_filename);
        EXPECT_EQ(result.out, std::string{});
        EXPECT_EQ(result.err, std::string{});
        RAPTOR_ASSERT_ZERO_EXIT(result);
    }

    ASSERT_TRUE(std::filesystem::exists(search_filename));
    search_result const result = parse_search_result(search_filename);

    EXPECT_EQ(result.user_bins.size(), user_bins.size());
    for (auto const & user_bin : result.user_bins)
    {
        std::string const start = test_directory.path().string() + "/user_bin_";
        ASSERT_TRUE(user_bin.starts_with(start)) << user_bin;
        EXPECT_LT(std::stoull(user_bin.substr(start.size())), user_bins.size()); // index is in correct range
    }

    EXPECT_EQ(result.queries.size(), queries.no_hits.size());
    for (size_t i = 0; i < result.queries.size(); ++i)
        EXPECT_TRUE(result.hits[i].empty()) << "Query: " << result.queries[i] << "is not empty";
}

INSTANTIATE_TEST_SUITE_P(
    search_partitioned_hibf_suite,
    search_partitioned_hibf,
    testing::Combine(testing::Values(2, 4, 8),
                     testing::Values(0, 1),
                     testing::Values("blocked", "sorted", "folded", "weighted_fold", "similarity")),
    [](testing::TestParamInfo<search_partitioned_hibf::ParamType> const & info)
    {
        std::string name = std::to_string(std::get<0>(info.param)) + "_parts" + std::to_string(std::get<1>(info.param))
                         + "_error" + std::get<2>(info.param) + "_approach";
        return name;
    });

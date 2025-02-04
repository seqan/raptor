// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <cassert>
#include <charconv>
#include <fstream>
#include <iostream>
#include <ranges>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include <sharg/parser.hpp>

#include <hibf/contrib/robin_hood.hpp>

struct parser_options
{
    std::filesystem::path result_file{};
    std::filesystem::path truth_file{};
    std::string result_format{};
    std::filesystem::path mantis_query_names{"none"};
    double mantis_threshold{-1.0};
};

void init_parser(sharg::parser & parser, parser_options & options)
{
    parser.add_option(options.result_file,
                      sharg::config{.short_id = '\0',
                                    .long_id = "results",
                                    .description = "The result file of a tool. Which tool is chosen by --result-format",
                                    .required = true,
                                    .validator = sharg::input_file_validator{}});
    parser.add_option(options.truth_file,
                      sharg::config{.short_id = '\0',
                                    .long_id = "truth",
                                    .description = "The ground truth to compare against.",
                                    .required = true,
                                    .validator = sharg::input_file_validator{}});
    parser.add_option(
        options.result_format,
        sharg::config{.short_id = '\0',
                      .long_id = "result-format",
                      .description = "The format of the result file.",
                      .required = true,
                      .validator = sharg::value_list_validator{"raptor", "mantis", "bifrost", "metagraph", "cobs"}});
    parser.add_option(options.mantis_query_names,
                      sharg::config{.short_id = '\0',
                                    .long_id = "mantis-query-names",
                                    .description = "If the format is mantis, the query ids are needed, in the same "
                                                   "order as the query file given to mantis when querying."});
    parser.add_option(
        options.mantis_threshold,
        sharg::config{.short_id = '\0',
                      .long_id = "mantis-threshold",
                      .description =
                          "If the format is mantis, the threshold given to all others tools when querying needed."});
}

auto create_truth_maps(std::filesystem::path const truth_file)
{
    // result data structures
    robin_hood::unordered_map<std::string, uint64_t> ub_name_to_id;
    robin_hood::unordered_map<std::string, std::vector<uint64_t>> hit_map;

    // parse header
    std::string line_buffer{};
    uint64_t buffer{};
    std::ifstream truth_file_in{truth_file};

    // Header contains lines: "some_number <tab> reference_name"
    std::cout << "[Truth] Parse header ..." << std::endl;
    while (std::getline(truth_file_in, line_buffer) && line_buffer.starts_with("#")
           && !line_buffer.starts_with("#QUERY_NAME"))
    {
        auto tab_it{line_buffer.begin() + line_buffer.find('\t')};
        std::string_view const id_value{line_buffer.begin() + 1 /* skip # */, tab_it};
        std::string_view const name_key{++tab_it, line_buffer.end()};
        std::from_chars(id_value.data(), id_value.data() + id_value.size(), buffer);
        ub_name_to_id.emplace(name_key, buffer);
    }

    assert(line_buffer == "#QUERY_NAME\tUSER_BINS");

    std::cout << "[Truth] Parse results ..." << std::endl;
    std::vector<uint64_t> result_user_bins;
    while (std::getline(truth_file_in, line_buffer))
    {
        result_user_bins.clear();

        auto tab_it{line_buffer.begin() + line_buffer.find('\t')};
        std::string const id{line_buffer.begin(), tab_it};
        std::string_view const bins{++tab_it, line_buffer.end()};

        for (auto && user_bin : bins | std::views::split(','))
        {
            // std::string_view id_value(user_bin.begin(), user_bin.end()); // doesn't work??
            std::string ub;
            for (char const chr : user_bin)
                ub.push_back(chr);

            std::from_chars(ub.data(), ub.data() + ub.size(), buffer);
            result_user_bins.push_back(buffer);
        }

        std::sort(result_user_bins.begin(), result_user_bins.end());

        hit_map[id] = result_user_bins;
    }

    std::cout << "Hit map contains " << hit_map.size() << " entries." << std::endl;

    return std::make_pair(ub_name_to_id, hit_map);
}

std::pair<size_t, size_t> compare_result_vectors(std::vector<uint64_t> const & truth,
                                                 std::vector<uint64_t> const & result,
                                                 size_t & FPS,
                                                 size_t & FNS)
{
    auto truth_hit = truth.begin();
    auto result_hit = result.begin();

    while (truth_hit != truth.end() && result_hit != result.end())
    {
        if (*truth_hit == *result_hit)
        {
            ++truth_hit;
            ++result_hit;
        }
        else
        {
            if (*truth_hit < *result_hit)
            {
                ++FNS;
                ++truth_hit;
            }
            else
            {
                ++FPS;
                ++result_hit;
            }
        }
    }
    while (truth_hit != truth.end())
    {
        ++FNS;
        ++truth_hit;
    }
    while (result_hit != result.end())
    {
        ++FPS;
        ++result_hit;
    }

    return {FPS, FNS};
}

void print_result(size_t const FPS, size_t const FNS)
{
    std::cout << std::endl;
    std::cout << "FPs " << FPS << std::endl;
    std::cout << "FNs " << FNS << std::endl;
    std::cout << std::endl;

    double const ALL = 25321.0 * 10012324.0;
    // ALL_P was obtained by counting hits in truth file
    // /project/archive-index-data/smehringer/raptor_bench/the_one_and_only.truth
    double const ALL_P = 855109966.0; // TP + FN
    double const ALL_N = ALL - ALL_P; // FP + TN

    double const TP = ALL_P - FNS;
    double const TN = ALL_N - FPS;
    double const PP = ALL_P - FNS + FPS;

    double const accuracy = ((TP + TN) * 100.0) / ALL;
    double const fpr = FPS / ALL_N;
    double const fdr = FPS / PP; // FPS / (All hits recorded as true by tool)
    double const sensitivity = TP / ALL_P;
    double const precision = TP / PP; // FPS / (All hits recorded as true by tool)
    double const specificity = TN / ALL_N;
    double const f1s = (2 * TP) / (2 * TP + FPS + FNS);

    std::cout << "Accuracy " << accuracy << std::endl;
    std::cout << "FPR " << fpr << std::endl;
    std::cout << "FDR " << fdr << std::endl;
    std::cout << "Sensitivity " << sensitivity << std::endl;
    std::cout << "Precision " << precision << std::endl;
    std::cout << "Specificity " << specificity << std::endl;
    std::cout << "F1 " << f1s << std::endl;
}

void print_progress(double const query_count, double const number_of_querries)
{
    int const barWidth = 70;
    double const progress = query_count / number_of_querries;

    std::cout << "[";
    int const pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i)
    {
        if (i < pos)
            std::cout << "=";
        else if (i == pos)
            std::cout << ">";
        else
            std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}

// =============================================================================
// RAPTOR
// =============================================================================

void compare_raptor_to_truth(std::filesystem::path const raptor_file_name,
                             robin_hood::unordered_map<std::string, uint64_t> const & truth_id_map,
                             robin_hood::unordered_map<std::string, std::vector<uint64_t>> const & hit_map)
{
    std::ifstream raptor_file_in{raptor_file_name};
    std::string raptor_line;
    size_t query_count{}; // in the end we will check if all queries are in the file

    size_t FPS = 0;
    size_t FNS = 0;

    // parse header and fill map that maps the "index used through out the raptor file" to "the official reference name"
    robin_hood::unordered_map<std::string, std::string> idx_to_id;

    std::cout << "[Raptor] Parse header ..." << std::endl;
    // Header first contains parameters starting with "##"
    while (std::getline(raptor_file_in, raptor_line) && raptor_line.starts_with("##"))
        ; // skip

    // Header contains lines: "some_number <tab> reference_name"
    do
    {
        auto tab_it{raptor_line.begin() + raptor_line.find('\t')};
        std::string_view const id_value{raptor_line.begin() + 1 /* skip # */, tab_it};
        std::string_view const name{++tab_it, raptor_line.end()};
        idx_to_id.emplace(id_value, name);
    }
    while (std::getline(raptor_file_in, raptor_line) && raptor_line.starts_with("#")
           && !raptor_line.starts_with("#QUERY_NAME"));

    assert(raptor_line == "#QUERY_NAME\tUSER_BINS");

    std::cout << "[Raptor] Parse results ..." << std::endl;
    std::vector<uint64_t> result_user_bins;
    while (std::getline(raptor_file_in, raptor_line))
    {
        // retrieve user bins
        result_user_bins.clear();

        auto tab_it{raptor_line.begin() + raptor_line.find('\t')};
        std::string id{raptor_line.begin(), tab_it};
        std::string_view const bins{++tab_it, raptor_line.end()};

        for (auto && user_bin : bins | std::views::split(','))
        {
            // std::string_view id_value(user_bin.begin(), user_bin.end()); // doesn't work??
            std::string ub;
            for (char const chr : user_bin)
                ub.push_back(chr);

            // -----------------------------------------------------------------
            // DEBUG OUTPUT
            // if (idx_to_id.find(ub) == idx_to_id.end())
            //     std::cout << "ub"<< ub << std::endl;
            // if (truth_id_map.find(idx_to_id.at(ub)) == truth_id_map.end())
            //     std::cout << idx_to_id.at(ub) << std::endl;
            // -----------------------------------------------------------------

            result_user_bins.push_back(truth_id_map.at(idx_to_id.at(ub)));
        }

        std::sort(result_user_bins.begin(), result_user_bins.end()); // compare script afterwards requires sorted UBs

        // compare vectors
        compare_result_vectors(hit_map.at(id), result_user_bins, FPS, FNS);

        ++query_count;
        print_progress(query_count, hit_map.size());
    }

    if (query_count != hit_map.size())
        throw std::runtime_error{"The result file did only contain " + std::to_string(query_count) + "/"
                                 + std::to_string(hit_map.size()) + " queries."};

    print_result(FPS, FNS);
}

// =============================================================================
// MANTIS
// =============================================================================

std::vector<std::string> parse_query_names(std::filesystem::path const & query_names_file)
{
    std::string line_buffer{};
    std::vector<std::string> query_names;
    std::ifstream query_names_in{query_names_file};

    std::cerr << "Reading " << query_names_file << " ... " << std::flush;
    // Contains lines: "query_name"
    while (std::getline(query_names_in, line_buffer))
        query_names.push_back(line_buffer);
    std::cerr << "Done" << std::endl;
    return query_names;
}

/* ## Threshold:
 * Let:
 *   * p patternsize
 *   * k k-mer size
 *   * e errors
 *   * c k-mer count
 *   * t threshold
 * Then:
 *   * c = p - k + 1 [Lemma A]
 *   * p = c + k - 1 [Eq 1]
 *   * t = p - (e + 1) * k + 1 [Lemma B]
 *   * t = c + k - 1 - (e + 1) * k + 1 [Eq 1 + Lemma B]
 *   * t = c + 1 - 1 + k - k - e * k
 *   * t = c - e * k
 * However, mantis counts unique kmers, not taking into account multiplicity.
 * Hence, we may further substract a constant from the threshold.
 */
class thresholder
{
private:
    // size_t const destroyed_kmers{};
    double percentage{};

public:
    thresholder() = default;
    thresholder(thresholder const &) = default;
    thresholder(thresholder &&) = default;
    thresholder & operator=(thresholder const &) = default;
    thresholder & operator=(thresholder &&) = default;
    ~thresholder() = default;

    explicit thresholder(double const perc) : percentage(perc)
    {}

    [[nodiscard]] constexpr size_t get(size_t const kmer_count) const noexcept
    {
        return static_cast<size_t>(percentage * kmer_count);
    }

    // explicit thresholder(config const & cfg) :
    //     destroyed_kmers(cfg.number_of_errors * cfg.kmer_size + cfg.threshold_grace)
    // {}

    // [[nodiscard]] constexpr size_t get(size_t const kmer_count) const noexcept
    // {
    //     return (kmer_count > destroyed_kmers) ? (kmer_count - destroyed_kmers) : 0u;
    // }
};

void compare_mantis_to_truth(std::filesystem::path const mantis_file_name,
                             std::filesystem::path const query_names_file,
                             double const the_given_threshold,
                             robin_hood::unordered_map<std::string, uint64_t> const & truth_id_map,
                             robin_hood::unordered_map<std::string, std::vector<uint64_t>> const & hit_map)
{
    std::ifstream mantis_file_in{mantis_file_name};
    std::string mantis_line;
    size_t query_count{};                             // in the end we will check if all queries are in the file
    thresholder const threshold{the_given_threshold}; // Helper for computing the threshold.
    size_t mantis_threshold{};                        // Needs to be set for each query.
    size_t current_query_number{};
    std::vector<uint64_t> results;
    std::string ub_name_buffer;

    size_t FPS = 0;
    size_t FNS = 0;

    // All query names
    std::vector<std::string> const query_names{parse_query_names(query_names_file)};

    // recreate a new truth_id_map becuase the ids path are from squeakr files not the original files
    robin_hood::unordered_map<std::string, uint64_t> truth_only_filename_id_map;

    for (auto && [key, value] : truth_id_map)
    {
        // truth ids: /path/v1/files/GCF_000019125.1_ASM1912v1_genomic.fna.gz
        std::string new_key{key.begin() + key.find_last_of('/') + 1, key.end() - 7};
        truth_only_filename_id_map[new_key] = value;
    }
    std::cout << truth_only_filename_id_map.begin()->first << std::endl;

    // ## Mantis results ##
    // There is no query ID, instead they are enumerated: seq0 - seqX
    // For each read:
    //   * seqX <tab> kmers in query
    //   * For each user bin that has kmers, the unique kmer count is reported:
    //     * absolute_path/to/ub.squeakr <tab> count
    // E.g.:
    // seq0    219
    // [...]/GCF_016028855.1_ASM1602885v1_genomic.squeakr       173
    // [...]/GCF_016128175.1_ASM1612817v1_genomic.squeakr       205
    // [...]/GCF_020162095.1_ASM2016209v1_genomic.squeakr       1
    // seq1    219
    auto parse_user_bin_id = [&ub_name_buffer, &truth_only_filename_id_map](std::string const & line)
    {
        ub_name_buffer.assign(line.begin() + line.find_last_of('/') + 1, line.begin() + line.find_last_of('.'));

        // -----------------------------------------------------------------
        // DEBUG OUTPUT
        if (truth_only_filename_id_map.find(ub_name_buffer) == truth_only_filename_id_map.end())
            std::cout << "ub_name_buffer" << ub_name_buffer << std::endl;
        // -----------------------------------------------------------------

        return truth_only_filename_id_map.at(ub_name_buffer);
    };

    auto parse_kmer_count = [](std::string const & line)
    {
        size_t result{};
        std::string_view const sv{line.begin() + line.find('\t') + 1 /* Skip seqX */, line.end()};
        std::from_chars(sv.data(), sv.data() + sv.size(), result);
        return result;
    };

    std::string current_query_name;

    auto process_results = [&]()
    {
        std::sort(results.begin(), results.end());

        // -----------------------------------------------------------------
        // DEBUG OUTPUT
        if (hit_map.find(current_query_name) == hit_map.end())
            std::cout << "current_query_name" << current_query_name << std::endl;
        // -----------------------------------------------------------------

        compare_result_vectors(hit_map.at(current_query_name), results, FPS, FNS);

        ++query_count;

        results.clear();
    };

    // First line.

    if (std::getline(mantis_file_in, mantis_line))
    {
        assert(mantis_line.starts_with("seq"));
        assert(current_query_number < query_names.size());
        current_query_name = query_names[current_query_number++];
        mantis_threshold = threshold.get(parse_kmer_count(mantis_line));
    }

    while (std::getline(mantis_file_in, mantis_line))
    {
        if (mantis_line.starts_with("seq")) // new query
        {
            // Process the results of the previous query.
            process_results();

            // Set new query name.
            assert(current_query_number < query_names.size());
            current_query_name = query_names[current_query_number++];

            // Compute threshold for current query.
            mantis_threshold = threshold.get(parse_kmer_count(mantis_line));
        }
        else
        {
            if (parse_kmer_count(mantis_line) >= mantis_threshold)
                results.push_back(parse_user_bin_id(mantis_line));
        }
    }

    // Write last results.
    process_results();

    if (query_count != hit_map.size())
        throw std::runtime_error{"The result file did only contain " + std::to_string(query_count) + "/"
                                 + std::to_string(hit_map.size()) + " queries."};

    print_result(FPS, FNS);
}

void compare_bifrost_to_truth(std::filesystem::path const bifrost_file_name,
                              robin_hood::unordered_map<std::string, uint64_t> const & truth_id_map,
                              robin_hood::unordered_map<std::string, std::vector<uint64_t>> const & hit_map)
{
    // Process bifrost results
    std::ifstream bifrost_result_in{bifrost_file_name};
    size_t query_count{};
    std::vector<uint64_t> results;
    size_t FPS = 0;
    size_t FNS = 0;

    // Buffers for file I/O
    std::string line_buffer{};

    std::vector<uint64_t> bifrost_user_bins{};

    // ## Bifrost results ##
    // Bifrost outputs a matrix. Column names = user bin id. Row names = query names
    auto split_line_by_tab_and = [](std::string_view bifrost_line, auto do_me)
    {
        std::string_view::size_type current_pos = 0;
        std::string_view::size_type tab_pos{bifrost_line.find('\t')};
        size_t column_idx{};

        while (tab_pos != std::string_view::npos)
        {
            auto current = std::string(&bifrost_line[current_pos], tab_pos - current_pos);
            do_me(current, column_idx);
            current_pos = tab_pos + 1;
            tab_pos = bifrost_line.find('\t', current_pos);
            ++column_idx;
        }
        // process last ub
        auto last = std::string(&bifrost_line[current_pos], bifrost_line.size() - current_pos);
        do_me(last, column_idx);
    };

    auto parse_header_user_bin_id = [&bifrost_user_bins, &truth_id_map](std::string const & sv, size_t idx)
    {
        if (idx != 0)
        {
            try
            {
                bifrost_user_bins.push_back(truth_id_map.at(sv));
            }
            catch (std::exception const & e)
            {
                std::cerr << "Could not find id: " << sv << std::endl;
                throw e;
            }
        }
        else
        {
            bifrost_user_bins.push_back(0); // don't mess up the indices
        }
    };

    auto insert_if_one = [&results, &bifrost_user_bins](std::string_view sv, size_t idx)
    {
        if (sv == std::string_view{"1"}) // excludes 0 and the first column which is alywas the query name
        {
            results.push_back(bifrost_user_bins[idx]);
        }
    };

    std::cerr << "Processing " << bifrost_file_name << " ... " << std::endl;

    // Parse header line
    if (std::getline(bifrost_result_in, line_buffer))
    {
        assert(line_buffer.starts_with("query_name"));
        split_line_by_tab_and(line_buffer, parse_header_user_bin_id);
    }

    assert(truth_id_map.size() + 1 == bifrost_user_bins.size());
    std::cerr << "Successfully parsed Header line ... " << std::endl;

    while (std::getline(bifrost_result_in, line_buffer))
    {
        results.clear();

        split_line_by_tab_and(line_buffer, insert_if_one);

        std::string_view::size_type tab_pos{line_buffer.find('\t')};
        std::string query_name(&line_buffer[0], tab_pos);

        compare_result_vectors(hit_map.at(query_name), results, FPS, FNS);

        ++query_count;
        print_progress(query_count, hit_map.size());
    }

    if (query_count != hit_map.size())
        throw std::runtime_error{"The result file did only contain " + std::to_string(query_count) + "/"
                                 + std::to_string(hit_map.size()) + " queries."};

    print_result(FPS, FNS);
}

void compare_metagraph_to_truth(std::filesystem::path const metagraph_file_name,
                                robin_hood::unordered_map<std::string, uint64_t> const & truth_id_map,
                                robin_hood::unordered_map<std::string, std::vector<uint64_t>> const & hit_map)
{

    // Process metagraph results
    std::ifstream metagraph_result_in{metagraph_file_name};
    size_t query_count{};
    std::vector<uint64_t> results;
    size_t FPS = 0;
    size_t FNS = 0;

    // Buffers for file I/O
    std::string metagraph_line{};

    while (std::getline(metagraph_result_in, metagraph_line))
    {
        results.clear();

        std::string_view::size_type first_tab{metagraph_line.find('\t')};
        std::string_view::size_type second_tab{metagraph_line.find('\t', first_tab + 1)};
        std::string query_name(metagraph_line.begin() + first_tab + 1, metagraph_line.begin() + second_tab);
        std::string bins{metagraph_line.begin() + second_tab + 1, metagraph_line.end()};

        for (auto && user_bin : bins | std::views::split(':'))
        {
            // std::string_view id_value(user_bin.begin(), user_bin.end()); // doesn't work??
            std::string ub;
            for (char const chr : user_bin)
                ub.push_back(chr);

            // -----------------------------------------------------------------
            // DEBUG OUTPUT
            // if (truth_id_map.find(ub) == truth_id_map.end())
            //     std::cout << ub << std::endl;
            // -----------------------------------------------------------------

            results.push_back(truth_id_map.at(ub));
        }

        // -----------------------------------------------------------------
        // DEBUG OUTPUT
        // if (hit_map.find(query_name) == hit_map.end())
        //     std::cout << query_name << std::endl;
        // -----------------------------------------------------------------

        compare_result_vectors(hit_map.at(query_name), results, FPS, FNS);

        ++query_count;
        print_progress(query_count, hit_map.size());
    }

    if (query_count != hit_map.size())
        throw std::runtime_error{"The result file did only contain " + std::to_string(query_count) + "/"
                                 + std::to_string(hit_map.size()) + " queries."};

    print_result(FPS, FNS);
}

void compare_cobs_to_truth(std::filesystem::path const cobs_file_name,
                           robin_hood::unordered_map<std::string, uint64_t> const & truth_id_map,
                           robin_hood::unordered_map<std::string, std::vector<uint64_t>> const & hit_map)
{
    // recreate a new truth_id_map becuase the ids are truncated in COBS output :(
    robin_hood::unordered_map<std::string, uint64_t> truth_trunc_id_map;

    for (auto && [key, value] : truth_id_map)
    {
        std::string new_key{key.begin() + key.find_last_of('/') + 1, key.begin() + key.find('.')};
        truth_trunc_id_map[new_key] = value;
    }

    std::ifstream cobs_result_in{cobs_file_name};
    size_t query_count{};
    std::vector<uint64_t> results;
    size_t FPS = 0;
    size_t FNS = 0;
    std::string cobs_line{};
    std::string current_query_name;

    auto parse_query_name = [&current_query_name](std::string_view const line)
    {
        current_query_name = std::string{line.begin() + 1 /* skip * */, line.begin() + line.find('\t')};
    };

    auto process_results = [&]()
    {
        std::sort(results.begin(), results.end());

        // -----------------------------------------------------------------
        // DEBUG OUTPUT
        // if (hit_map.find(current_query_name) == hit_map.end())
        //     std::cout << "current_query_name: " << current_query_name << std::endl;
        // -----------------------------------------------------------------

        compare_result_vectors(hit_map.at(current_query_name), results, FPS, FNS);

        results.clear();
    };

    auto parse_user_bin_id = [&truth_trunc_id_map](std::string_view const line)
    {
        std::string ub{line.begin(), line.begin() + line.find('\t')};

        // -----------------------------------------------------------------
        // DEBUG OUTPUT
        // if (truth_trunc_id_map.find(ub) == truth_trunc_id_map.end())
        //     std::cout << "UB:" << ub << std::endl;
        // -----------------------------------------------------------------

        return truth_trunc_id_map.at(ub);
    };

    // First line.
    /* COBS result looks like this:
     * ```
     * *query_name1  [some number]
     * reference_name1   [COUNT]
     * reference_name2   [COUNT]
     * ...
     * ```
     */
    if (std::getline(cobs_result_in, cobs_line))
    {
        assert(cobs_line.starts_with("*"));
        parse_query_name(cobs_line);
    }

    while (std::getline(cobs_result_in, cobs_line))
    {
        if (cobs_line.starts_with("*")) // new query
        {
            // Process the results of the previous query.
            process_results();

            // Set new query name.
            parse_query_name(cobs_line);

            ++query_count;
            print_progress(query_count, hit_map.size());
        }
        else
        {
            results.push_back(parse_user_bin_id(cobs_line));
        }
    }

    // Write last results.
    process_results();
    ++query_count;

    if (query_count != hit_map.size())
        throw std::runtime_error{"The result file did only contain " + std::to_string(query_count) + "/"
                                 + std::to_string(hit_map.size()) + " queries."};

    print_result(FPS, FNS);
}

int main(int argc, char ** argv)
{
    std::cout << "Hello..." << std::endl;

    sharg::parser parser{"compare_to_truth", argc, argv, sharg::update_notifications::off};
    parser.info.author = "Svenja Mehringer, Enrico Seiler";
    parser.info.email = "enrico.seiler@fu-berlin.de";
    parser.info.short_description = "Compares mantis and raptor results.";
    parser.info.version = "0.0.1";

    parser_options options{};
    init_parser(parser, options);

    try
    {
        parser.parse();

        if (options.result_format == "mantis" && options.mantis_query_names == "none")
            throw sharg::parser_error{"For mantis results you need to pass the query ids."};

        if (options.result_format == "mantis" && options.mantis_threshold == -1.0)
            throw sharg::parser_error{"For mantis results you need to pass the threshold."};
    }
    catch (sharg::parser_error const & ext)
    {
        std::cerr << "[Error] " << ext.what() << '\n';
        std::exit(-1);
    }

    std::cout << "Create truth data structures to compare against..." << std::endl;
    auto [id_map, hit_map] = create_truth_maps(options.truth_file);

    std::cout << "Compare " << options.result_format << " result..." << std::endl;

    if (options.result_format == "raptor")
        compare_raptor_to_truth(options.result_file, id_map, hit_map);
    else if (options.result_format == "mantis")
        compare_mantis_to_truth(options.result_file,
                                options.mantis_query_names,
                                options.mantis_threshold,
                                id_map,
                                hit_map);
    else if (options.result_format == "bifrost")
        compare_bifrost_to_truth(options.result_file, id_map, hit_map);
    else if (options.result_format == "metagraph")
        compare_metagraph_to_truth(options.result_file, id_map, hit_map);
    else if (options.result_format == "cobs")
        compare_cobs_to_truth(options.result_file, id_map, hit_map);
}

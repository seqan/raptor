// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/eseiler/minimizer_thresholds/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides stuff.
 */

#include <chrono>
#include <seqan3/std/filesystem>
#include <numeric>
#include <unordered_set>
#include <set>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/views/async_input_buffer.hpp>

#include <cereal/types/vector.hpp>

#include <heuristic_threshold.hpp>
#include <lemma_threshold.hpp>
#include <minimizer_model.hpp>

struct my_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

struct cmd_arguments
{
    std::filesystem::path reference_file{};
    std::filesystem::path query_file{};
    std::filesystem::path output_directory{};
    uint64_t window_size{26};
    uint8_t kmer_size{20};
    uint8_t errors{3};
    uint64_t pattern_size{};
    double tau{0.99};
    bool from_file{};
    std::vector<std::string> method{};
};

struct threshold_result
{
    std::vector<size_t> threshold_per_read;
    size_t number_of_hits{};
    std::chrono::duration<double, std::milli> precompute_time{0.0};
    std::chrono::duration<double, std::milli> threshold_time{0.0};
    std::string method;
};

void do_cerealisation_out(std::vector<size_t> const & vec, cmd_arguments const & args, threshold_result const & result)
{
    std::filesystem::path filename = args.output_directory / ("binary_" + result.method +
                                                              "_w" + std::to_string(args.window_size) +
                                                              "_k" + std::to_string(args.kmer_size) +
                                                              "_e" + std::to_string(args.errors) +
                                                              "_tau" + std::to_string(args.tau));
    std::ofstream os{filename, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(vec);
}

void do_cerealisation_in(std::vector<size_t> & vec, cmd_arguments const & args, threshold_result const & result)
{
    std::filesystem::path filename = args.output_directory / ("binary_" + result.method +
                                                              "_w" + std::to_string(args.window_size) +
                                                              "_k" + std::to_string(args.kmer_size) +
                                                              "_e" + std::to_string(args.errors) +
                                                              "_tau" + std::to_string(args.tau));
    std::ifstream is{filename, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    iarchive(vec);
}

std::unordered_set<uint64_t> hash_reference(minimizer & mini, std::filesystem::path const & reference_file)
{
    std::unordered_set<uint64_t> reference_hashes{};
    seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> reference_in{reference_file};
    for (auto & [seq] : reference_in)
    {
        mini.compute(seq);
        reference_hashes.insert(mini.minimizer_hash.begin(), mini.minimizer_hash.end());
    }
    return reference_hashes;
}

void print_results(threshold_result const & result)
{
    if (std::ranges::empty(result.threshold_per_read))
    {
        seqan3::debug_stream << "Minimal threshold: 0\n";
        seqan3::debug_stream << "Average threshold: 0\n";
    }
    else
    {
        seqan3::debug_stream << "Minimal threshold: "
                             << *std::min_element(result.threshold_per_read.begin(), result.threshold_per_read.end())
                             << '\n';
        seqan3::debug_stream << "Average threshold: "
                             << static_cast<double>(std::accumulate(result.threshold_per_read.begin(),
                                                                    result.threshold_per_read.end(), 0))
                                / result.threshold_per_read.size()
                            << '\n';
    }
    seqan3::debug_stream << "Hits: " << result.number_of_hits << '\n';
    seqan3::debug_stream << "Precompute time: " << result.precompute_time.count() << '\n';
    seqan3::debug_stream << "Threshold time: " << result.threshold_time.count() << '\n';
    // seqan3::debug_stream << '\n';
}

void write_results(cmd_arguments const & args,
                   threshold_result const & result)
{
    std::filesystem::path out_file = args.output_directory / ("result_" + result.method +
                                                              "_w" + std::to_string(args.window_size) +
                                                              "_k" + std::to_string(args.kmer_size) +
                                                              "_e" + std::to_string(args.errors) +
                                                              "_tau" + std::to_string(args.tau) +
                                                              ".csv");

    std::ofstream out{out_file};
    if (!out.is_open())
        throw std::filesystem::filesystem_error("Could not open file for writing.",
                                                out_file,
                                                std::make_error_code(std::errc::bad_file_descriptor));

    out << "method,min_threshold,mean_threshold,hits,precompute_time,threshold_time\n";
    out << result.method << ',';

    if (std::ranges::empty(result.threshold_per_read))
    {
        out << "0,0,";
    }
    else
    {
        out << *std::min_element(result.threshold_per_read.begin(), result.threshold_per_read.end())
            << ','
            << (static_cast<double>(std::accumulate(result.threshold_per_read.begin(),
                                    result.threshold_per_read.end(), 0))
               / result.threshold_per_read.size())
            << ',';
    }

    out << result.number_of_hits << ',';
    out << result.precompute_time.count() << ',';
    out << result.threshold_time.count() << '\n';
}

void compute_heuristic_threshold(cmd_arguments const & args)
{
    minimizer mini{window{args.window_size}, kmer{args.kmer_size}};
    std::unordered_set<uint64_t> reference_hashes = hash_reference(mini, args.reference_file);
    threshold_result result;
    result.method = "heuristic";
    std::vector<uint64_t> threshold_per_count(args.pattern_size, 0);

    // Process the queries.
    seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> query_in{args.query_file};
    for (auto & [seq] : query_in | seqan3::views::async_input_buffer(16))
    {
        size_t minimizer_count = 0;

        for (auto && hash : mini.hashes(seq))
            minimizer_count += reference_hashes.count(hash);

        uint64_t count = mini.minimizer_hash.size();

        auto start = std::chrono::high_resolution_clock::now();
        heuristic_threshold heuristic(mini);
        auto threshold = heuristic.threshold(args.errors);
        auto end = std::chrono::high_resolution_clock::now();
        threshold_per_count[count] = threshold;
        result.threshold_time += end - start;
        result.threshold_per_read.push_back(threshold);
        result.number_of_hits += (threshold <= minimizer_count);
    }

    print_results(result);
    write_results(args, result);
    seqan3::debug_stream << threshold_per_count << '\n';
    seqan3::debug_stream << '\n';
}

void compute_lemma_threshold(cmd_arguments const & args)
{
    minimizer mini{window{args.window_size}, kmer{args.kmer_size}};
    std::unordered_set<uint64_t> reference_hashes = hash_reference(mini, args.reference_file);
    threshold_result result;
    result.method = "lemma";

    // Process the queries.
    seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> query_in{args.query_file};
    for (auto & [seq] : query_in | seqan3::views::async_input_buffer(16))
    {
        size_t minimizer_count = 0;

        for (auto && hash : mini.hashes_multi(seq))
            minimizer_count += reference_hashes.count(hash);

        auto start = std::chrono::high_resolution_clock::now();
        lemma_threshold lemma{};
        auto threshold = lemma.threshold(args.pattern_size, args.window_size, args.errors);
        auto end = std::chrono::high_resolution_clock::now();
        result.threshold_time += end - start;
        result.threshold_per_read.push_back(threshold);
        result.number_of_hits += (threshold <= minimizer_count);
    }

    print_results(result);
    write_results(args, result);
}

void compute_simple_model(cmd_arguments const & args,
                          bool const indirect,
                          bool const overlapping)
{
    minimizer mini{window{args.window_size}, kmer{args.kmer_size}};
    std::unordered_set<uint64_t> reference_hashes = hash_reference(mini, args.reference_file);
    threshold_result result;
    result.method = (indirect && overlapping ? "indirect_overlapping" :
                    (indirect ? "indirect" :
                    (overlapping ? "overlapping" :
                    "simple")));

    std::vector<uint64_t> threshold_per_count(args.pattern_size, 0);
    size_t kmers_per_window = args.window_size - args.kmer_size + 1;
    size_t kmers_per_pattern = args.pattern_size - args.kmer_size + 1;
    size_t minimal_number_of_minimizers = kmers_per_window == 1 ? kmers_per_pattern :
                                                                  kmers_per_pattern / (kmers_per_window - 1);
    size_t maximal_number_of_minimizers = args.pattern_size - args.window_size + 1;

    std::vector<size_t> precomp_thresholds;
    auto start2 = std::chrono::high_resolution_clock::now();

    if (args.from_file)
    {
        do_cerealisation_in(precomp_thresholds, args, result);
    }
    else
    {
        precomp_thresholds = precompute_threshold(args.pattern_size,
                                                  args.window_size,
                                                  args.kmer_size,
                                                  args.errors,
                                                  args.tau,
                                                  indirect,
                                                  overlapping);
    }

    auto end2 = std::chrono::high_resolution_clock::now();
    result.precompute_time += end2 - start2;
    if (!args.from_file)
        do_cerealisation_out(precomp_thresholds, args, result);

    // Process the queries.
    seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> query_in{args.query_file};
    for (auto & [seq] : query_in | seqan3::views::async_input_buffer(16))
    {
        size_t minimizer_count = 0;

        for (auto && hash : mini.hashes(seq))
            minimizer_count += reference_hashes.count(hash);

        uint64_t count = mini.minimizer_hash.size();

        auto start = std::chrono::high_resolution_clock::now();
        size_t index = std::min(mini.minimizer_hash.size() < minimal_number_of_minimizers ? 0 :
                                    mini.minimizer_hash.size() - minimal_number_of_minimizers,
                                maximal_number_of_minimizers - minimal_number_of_minimizers);
        auto threshold = precomp_thresholds[index];
        auto end = std::chrono::high_resolution_clock::now();
        threshold_per_count[count] = threshold;
        result.threshold_time += end - start;
        result.threshold_per_read.push_back(threshold);
        result.number_of_hits += (threshold <= minimizer_count);
    }

    print_results(result);
    write_results(args, result);
    seqan3::debug_stream << threshold_per_count << '\n';
    seqan3::debug_stream << '\n';
}

void initialize_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "Enrico Seiler";
    parser.info.short_description = "This programs computes minimizer thresholds using different approaches.";
    parser.info.version = "1.0.0";

    parser.add_option(args.reference_file, 'r', "reference", "Provide a reference file.", seqan3::option_spec::required,
                      seqan3::input_file_validator<seqan3::sequence_file_input<>>{});
    parser.add_option(args.query_file, 'q', "query", "Provide a query file.", seqan3::option_spec::required,
                      seqan3::input_file_validator<seqan3::sequence_file_input<>>{});
    parser.add_option(args.output_directory, 'o', "out", "", seqan3::option_spec::required,
                      seqan3::output_directory_validator{});
    parser.add_option(args.window_size, 'w', "window", "Choose the window size.", seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{1, 1000});
    parser.add_option(args.kmer_size, 'k', "kmer", "Choose the kmer size.", seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{1, 32});
    parser.add_option(args.errors, 'e', "error", "Choose the number of errors.", seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{0, 5});
    parser.add_option(args.tau, 't', "tau", "Threshold for probabilistic models.", seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{0, 1});
    parser.add_option(args.pattern_size, 'p', "pattern", "Choose the pattern size. Only needed for methods other than "
                      "heuristic or lemma. Default: Use median of sequence lengths in query file.");
    parser.add_option(args.method, 'm', "method", "Choose the methods to compute the trheshold.",
                      seqan3::option_spec::required, seqan3::value_list_validator{"heuristic", "lemma", "simple",
                      "indirect", "overlap", "indirect-overlap", "all"});
    parser.add_flag(args.from_file, 'f', "from_file", "Load precomputed threshold from disk. Program must have "
                      "been run without this flag before.");
}

int main(int argc, char ** argv)
{
    seqan3::argument_parser myparser{"Minimizers", argc, argv, seqan3::update_notifications::off};
    cmd_arguments args{};

    initialize_argument_parser(myparser, args);

    try
    {
         myparser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        std::cerr << "[Error] " << ext.what() << "\n";
        return -1;
    }

    if (!args.pattern_size)
    {
        std::vector<uint64_t> sequence_lengths{};
        seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> query_in{args.query_file};
        for (auto & [seq] : query_in | seqan3::views::async_input_buffer(16))
        {
            sequence_lengths.push_back(std::ranges::size(seq));
        }
        std::sort(sequence_lengths.begin(), sequence_lengths.end());
        args.pattern_size = sequence_lengths[sequence_lengths.size()/2];
    }

    bool all = std::find(args.method.begin(), args.method.end(), "all") != args.method.end();
    bool heuristic = std::find(args.method.begin(), args.method.end(), "heuristic") != args.method.end();
    bool lemma = std::find(args.method.begin(), args.method.end(), "lemma") != args.method.end();
    bool simple = std::find(args.method.begin(), args.method.end(), "simple") != args.method.end();
    bool indirect = std::find(args.method.begin(), args.method.end(), "indirect") != args.method.end();
    bool overlap = std::find(args.method.begin(), args.method.end(), "overlap") != args.method.end();
    bool indirect_overlap = std::find(args.method.begin(), args.method.end(), "indirect-overlap") != args.method.end();

    if (all || heuristic)
    {
        seqan3::debug_stream << "Heuristic\n";
        compute_heuristic_threshold(args);
    }
    if (all || lemma)
    {
        seqan3::debug_stream << "Lemma\n";
        compute_lemma_threshold(args);
    }
    if (all || simple)
    {
        seqan3::debug_stream << "Model\n";
        compute_simple_model(args, false, false);
    }
    if (all || indirect)
    {
        seqan3::debug_stream << "Model with indirect\n";
        compute_simple_model(args, true, false);
    }
    if (all || overlap)
    {
        seqan3::debug_stream << "Model with overlapping\n";
        compute_simple_model(args, false, true);
    }
    if (all || indirect_overlap)
    {
        seqan3::debug_stream << "Model with indirect and overlapping\n";
        compute_simple_model(args, true, true);
    }

    return 0;
}

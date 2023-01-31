// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <filesystem>
#include <mutex>
#include <vector>

#include <robin_hood.h>

#include <sharg/all.hpp>

#include <seqan3/core/algorithm/detail/execution_handler_parallel.hpp>
#include <seqan3/io/views/async_input_buffer.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/utility/views/chunk.hpp>

#include <raptor/adjust_seed.hpp>
#include <raptor/argument_parsing/search_arguments.hpp>
#include <raptor/argument_parsing/validators.hpp>
#include <raptor/dna4_traits.hpp>
#include <raptor/search/do_parallel.hpp>
#include <raptor/threshold/precompute_correction.hpp>
#include <raptor/threshold/precompute_threshold.hpp>

void threshold_info(raptor::search_arguments const & arguments, std::string const & shape_string)
{
    double compute_time{};

    uint8_t const kmer_size{arguments.shape.size()};
    size_t const kmers_per_window = arguments.window_size - kmer_size + 1;
    size_t const kmers_per_pattern = arguments.query_length - kmer_size + 1;
    size_t const minimal_number_of_minimizers = kmers_per_pattern / kmers_per_window;
    size_t const maximal_number_of_minimizers = arguments.query_length - arguments.window_size + 1;

    auto const parameters = arguments.make_threshold_parameters();
    std::vector<size_t> const precomp_correction = raptor::threshold::precompute_correction(parameters);
    std::vector<size_t> const precomp_thresholds = raptor::threshold::precompute_threshold(parameters);

    seqan3::sequence_file_input<raptor::dna4_traits, seqan3::fields<seqan3::field::seq>> fin{arguments.query_file};
    using record_type = typename decltype(fin)::record_type;
    std::vector<record_type> records{};

    std::vector<size_t> minimiser_frequencies(maximal_number_of_minimizers + 1);
    std::mutex min_freq_mutex;

    auto worker = [&](size_t const start, size_t const end)
    {
        auto hash_adaptor = seqan3::views::minimiser_hash(seqan3::ungapped{kmer_size},
                                                          seqan3::window_size{arguments.window_size},
                                                          seqan3::seed{raptor::adjust_seed(kmer_size)});

        for (auto && [seq] : records | seqan3::views::slice(start, end))
        {
            size_t const minimiser_count{static_cast<size_t>(std::ranges::distance(seq | hash_adaptor))};

            {
                std::lock_guard const lock{min_freq_mutex};
                ++minimiser_frequencies[minimiser_count];
            }
        }
    };

    for (auto && record_batch : fin | seqan3::views::chunk((1ULL << 20)))
    {
        records.clear();
        std::ranges::move(record_batch, std::back_inserter(records));

        raptor::do_parallel(worker, records.size(), arguments.threads, compute_time);
    }

    std::ofstream out{arguments.out_file};

    // clang-format off
    out << "#query: " << arguments.query_file << '\n'
        << "#output: " << arguments.out_file << '\n'
        << "#kmer: " << std::to_string(kmer_size) << '\n'
        << "#shape: " << shape_string << '\n'
        << "#window: " << arguments.window_size << '\n'
        << "#error: " << std::to_string(arguments.errors) << '\n'
        << "#tau: " << arguments.tau << '\n'
        << "#p_max: " << arguments.p_max << '\n'
        << "#fpr: " << arguments.fpr << '\n'
        << "#pattern: " << arguments.query_length << '\n'
        << "#threads: " << std::to_string(arguments.threads) << '\n'
        << "##minimal_number_of_minimizers: " << minimal_number_of_minimizers << '\n'
        << "##maximal_number_of_minimizers: " << maximal_number_of_minimizers << '\n'
        << "##x: Number of minimizers\n"
        << "###x: Number reads with x minimizers\n"
        << "##t(x): Total threshold = t_p(x) + t_c(x)\n"
        << "##t_p(x): Probabilistic threshold\n"
        << "##t_c(x): Correction term\n"
        << 'x' << ','
        << "#x" << ','
        << "t(x)" << ','
        << "t_p(x)" << ','
        << "t_c(x)" << '\n';
    // clang-format on

    if (minimiser_frequencies.empty())
        return;

    size_t last_non_zero_index{minimiser_frequencies.size() - 1};
    while (last_non_zero_index && !minimiser_frequencies[last_non_zero_index])
        --last_non_zero_index;
    minimiser_frequencies.resize(last_non_zero_index + 1);

    for (size_t i{minimal_number_of_minimizers}; i <= last_non_zero_index; ++i)
    {
        size_t const minimiser_count{minimiser_frequencies[i]};
        if (!minimiser_count)
            continue;
        size_t const index{i - minimal_number_of_minimizers};
        // clang-format off
        out << i << ','
            << minimiser_count << ','
            << precomp_thresholds[index] + precomp_correction[index] << ','
            << precomp_thresholds[index] << ','
            << precomp_correction[index] << '\n';
        // clang-format on
    }
}

void init_search_parser(sharg::parser & parser, raptor::search_arguments & arguments, std::string & shape_string)
{
    arguments.cache_thresholds = false;

    parser.add_option(arguments.query_file,
                      sharg::config{.short_id = '\0',
                                    .long_id = "query",
                                    .description = "Provide a path to the query file.",
                                    .required = true,
                                    .validator = sharg::input_file_validator{}});
    parser.add_option(arguments.out_file,
                      sharg::config{.short_id = '\0',
                                    .long_id = "output",
                                    .description = "Provide a path to the output.",
                                    .required = true});
    parser.add_option(arguments.errors,
                      sharg::config{.short_id = '\0',
                                    .long_id = "error",
                                    .description = "The number of errors",
                                    .validator = raptor::positive_integer_validator{true}});
    parser.add_option(arguments.shape_size,
                      sharg::config{.short_id = '\0',
                                    .long_id = "kmer",
                                    .description = "The k-mer size. Mutually exclusive with --shape.",
                                    .validator = sharg::arithmetic_range_validator{1, 32}});
    parser.add_option(arguments.window_size,
                      sharg::config{.short_id = '\0',
                                    .long_id = "window",
                                    .description = "The window size.",
                                    .validator = raptor::positive_integer_validator{}});
    parser.add_option(shape_string,
                      sharg::config{.short_id = '\0',
                                    .long_id = "shape",
                                    .description = "The shape to use for k-mers. Mutually exclusive with --kmer.",
                                    .validator = sharg::regex_validator{"[01]+"}});
    parser.add_option(
        arguments.tau,
        sharg::config{.short_id = '\0',
                      .long_id = "tau",
                      .description = "Used in the dynamic thresholding. The higher tau, the lower the threshold.",
                      .validator = sharg::arithmetic_range_validator{0, 1}});
    parser.add_option(
        arguments.threshold,
        sharg::config{.short_id = '\0',
                      .long_id = "threshold",
                      .description = "If set, this threshold is used instead of the probabilistic models.",
                      .validator = sharg::arithmetic_range_validator{0, 1}});
    parser.add_option(
        arguments.p_max,
        sharg::config{.short_id = '\0',
                      .long_id = "p_max",
                      .description = "Used in the dynamic thresholding. The higher p_max, the lower the threshold.",
                      .validator = sharg::arithmetic_range_validator{0, 1}});
    parser.add_option(arguments.fpr,
                      sharg::config{.short_id = '\0',
                                    .long_id = "fpr",
                                    .description = "The false positive rate used for building the index.",
                                    .validator = sharg::arithmetic_range_validator{0, 1}});
    parser.add_option(arguments.query_length,
                      sharg::config{
                          .short_id = '\0',
                          .long_id = "query_length",
                          .description = "The pattern size. Default: Use median of sequence lengths in query file.",
                      });
    parser.add_option(arguments.threads,
                      sharg::config{.short_id = '\0',
                                    .long_id = "threads",
                                    .description = "The number of threads to use.",
                                    .validator = raptor::positive_integer_validator{}});
}

int main(int argc, char ** argv)
{
    sharg::parser parser{"threshold_info", argc, argv, sharg::update_notifications::off};
    parser.info.author = "Enrico Seiler";
    parser.info.author = "enrico.seiler@fu-berlin.de";
    parser.info.short_description = "Print thresholds.";
    parser.info.version = "0.0.1";

    raptor::search_arguments arguments{};
    std::string shape_string{};
    init_search_parser(parser, arguments, shape_string);

    try
    {
        parser.parse();
    }
    catch (sharg::parser_error const & ext)
    {
        std::cerr << "[Error] " << ext.what() << '\n';
        std::exit(-1);
    }

    if (parser.is_option_set("shape"))
    {
        if (parser.is_option_set("kmer"))
            throw sharg::parser_error{"You cannot set both shape and k-mer arguments."};

        uint64_t tmp{};

        std::from_chars(shape_string.data(), shape_string.data() + shape_string.size(), tmp, 2);

        arguments.shape = seqan3::shape{seqan3::bin_literal{tmp}};
    }
    else
    {
        arguments.shape = seqan3::shape{seqan3::ungapped{arguments.shape_size}};
        shape_string.assign(arguments.shape_size, '1');
    }

    std::filesystem::path output_directory = arguments.out_file.parent_path();
    std::error_code ec{};
    std::filesystem::create_directories(output_directory, ec);

    if (!output_directory.empty() && ec)
        throw sharg::parser_error{
            seqan3::detail::to_string("Failed to create directory\"", output_directory.c_str(), "\": ", ec.message())};

    if (!arguments.query_length)
    {
        std::vector<uint64_t> sequence_lengths{};
        seqan3::sequence_file_input<raptor::dna4_traits, seqan3::fields<seqan3::field::seq>> fin{arguments.query_file};
        for (auto & [seq] : fin | seqan3::views::async_input_buffer(16))
        {
            sequence_lengths.push_back(std::ranges::size(seq));
        }
        std::sort(sequence_lengths.begin(), sequence_lengths.end());
        arguments.query_length = sequence_lengths[sequence_lengths.size() / 2];
    }

    threshold_info(arguments, shape_string);
}

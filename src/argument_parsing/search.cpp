// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <seqan3/io/views/async_input_buffer.hpp>

#include <raptor/argument_parsing/search.hpp>
#include <raptor/search/search.hpp>

namespace raptor
{

void init_search_parser(seqan3::argument_parser & parser, search_arguments & arguments)
{
    init_shared_meta(parser);
    init_shared_options(parser, arguments);
    parser.add_option(arguments.ibf_file,
                      '\0',
                      "index",
                      arguments.is_socks ? "Provide a valid path to an IBF." :
                                           "Provide a valid path to an IBF. Parts: Without suffix _0",
                      seqan3::option_spec::required);
    parser.add_option(arguments.query_file,
                      '\0',
                      "query",
                      "Provide a path to the query file.",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{});
    parser.add_option(arguments.out_file,
                      '\0',
                      "output",
                      "Please provide a valid path to the output.",
                      seqan3::option_spec::required,
                      seqan3::output_file_validator{});
    parser.add_option(arguments.errors,
                      '\0',
                      "error",
                      "Choose the number of errors.",
                      arguments.is_socks ? seqan3::option_spec::hidden : seqan3::option_spec::standard,
                      positive_integer_validator{true});
    parser.add_option(arguments.tau,
                      '\0',
                      "tau",
                      "Threshold for probabilistic models.",
                      arguments.is_socks ? seqan3::option_spec::hidden : seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{0, 1});
    parser.add_option(arguments.threshold,
                      '\0',
                      "threshold",
                      "If set, this threshold is used instead of the probabilistic models.",
                      arguments.is_socks ? seqan3::option_spec::hidden : seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{0, 1});
    parser.add_option(arguments.pattern_size,
                      '\0',
                      "pattern",
                      "Choose the pattern size. Default: Use median of sequence lengths in query file.",
                      arguments.is_socks ? seqan3::option_spec::hidden : seqan3::option_spec::standard);
    parser.add_flag(arguments.compressed,
                    '\0',
                    "compressed",
                    "Build a compressed IBF.");
    parser.add_flag(arguments.write_time,
                    '\0',
                    "time",
                    "Write timing file.",
                    seqan3::option_spec::advanced);
}

void run_search(seqan3::argument_parser & parser, bool const is_socks)
{
    search_arguments arguments{};
    arguments.is_socks = is_socks;
    init_search_parser(parser, arguments);
    try_parsing(parser);

    // ==========================================
    // Various checks.
    // ==========================================

    if (!arguments.is_socks)
    {
        seqan3::input_file_validator<seqan3::sequence_file_input<>>{}(arguments.query_file);
    }

    arguments.treshold_was_set = parser.is_option_set("threshold");

    if (arguments.parts == 1)
    {
        seqan3::input_file_validator{}(arguments.ibf_file);
    }
    else
    {
        seqan3::input_file_validator validator{};
        for (size_t part{0}; part < arguments.parts; ++part)
        {
            validator(arguments.ibf_file.string() + std::string{"_"} + std::to_string(part));
        }
    }

    // ==========================================
    // Process --pattern.
    // ==========================================
    if (!arguments.is_socks && !arguments.pattern_size)
    {
        std::vector<uint64_t> sequence_lengths{};
        seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>> query_in{arguments.query_file};
        for (auto & [seq] : query_in | seqan3::views::async_input_buffer(16))
        {
            sequence_lengths.push_back(std::ranges::size(seq));
        }
        std::sort(sequence_lengths.begin(), sequence_lengths.end());
        arguments.pattern_size = sequence_lengths[sequence_lengths.size()/2];
    }

    // ==========================================
    // Read window and kmer size, and the bin paths.
    // ==========================================
    {
        std::ifstream is{arguments.parts == 1 ? arguments.ibf_file.string() :
                                                arguments.ibf_file.string() + std::string{"_0"},
                         std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        iarchive(arguments.kmer_size);
        iarchive(arguments.window_size);
        iarchive(arguments.bin_path);
        if (arguments.is_socks)
            arguments.pattern_size = arguments.kmer_size;
    }

    // ==========================================
    // Dispatch
    // ==========================================
    raptor_search(arguments);
};

} // namespace raptor

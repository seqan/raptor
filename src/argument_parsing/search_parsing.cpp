// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <seqan3/io/views/async_input_buffer.hpp>

#include <raptor/argument_parsing/init_shared_meta.hpp>
#include <raptor/argument_parsing/search_parsing.hpp>
#include <raptor/argument_parsing/validators.hpp>
#include <raptor/dna4_traits.hpp>
#include <raptor/index.hpp>
#include <raptor/search/search.hpp>

namespace raptor
{

void init_search_parser(seqan3::argument_parser & parser, search_arguments & arguments)
{
    init_shared_meta(parser);
    parser.info.examples = {"raptor search --error 2 --index raptor.index --query queries.fastq --output search.output"};
    parser.add_option(arguments.index_file,
                      '\0',
                      "index",
                      arguments.is_socks ? "Provide a valid path to an index." :
                                           "Provide a valid path to an index. Parts: Without suffix _0",
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
                      "Provide a path to the output.",
                      seqan3::option_spec::required);
    parser.add_option(arguments.errors,
                      '\0',
                      "error",
                      "The number of errors",
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
    parser.add_option(arguments.p_max,
                      '\0',
                      "p_max",
                      "Correction.",
                      arguments.is_socks ? seqan3::option_spec::hidden : seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{0, 1});
    parser.add_option(arguments.fpr,
                      '\0',
                      "fpr",
                      "fpr.",
                      arguments.is_socks ? seqan3::option_spec::hidden : seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{0, 1});
    parser.add_option(arguments.pattern_size,
                      '\0',
                      "pattern",
                      "The pattern size. Default: Use median of sequence lengths in query file.",
                      arguments.is_socks ? seqan3::option_spec::hidden : seqan3::option_spec::standard);
    parser.add_option(arguments.threads,
                      '\0',
                      "threads",
                      "The numer of threads to use.",
                      seqan3::option_spec::standard,
                      positive_integer_validator{});
    parser.add_flag(arguments.is_hibf,
                    '\0',
                    "hibf",
                    "Index is an HIBF.",
                    seqan3::option_spec::advanced);
    parser.add_flag(arguments.write_time,
                    '\0',
                    "time",
                    "Write timing file.",
                    seqan3::option_spec::advanced);
}

void search_parsing(seqan3::argument_parser & parser, bool const is_socks)
{
    search_arguments arguments{};
    arguments.is_socks = is_socks;
    init_search_parser(parser, arguments);
    parser.parse();

    // ==========================================
    // Various checks.
    // ==========================================

    std::filesystem::path output_directory = arguments.out_file.parent_path();
    std::error_code ec{};
    std::filesystem::create_directories(output_directory, ec);

// LCOV_EXCL_START
    if (!output_directory.empty() && ec)
        throw seqan3::argument_parser_error{seqan3::detail::to_string("Failed to create directory\"",
                                                                      output_directory.c_str(),
                                                                      "\": ",
                                                                      ec.message())};
// LCOV_EXCL_STOP

    if (!arguments.is_socks)
    {
        seqan3::input_file_validator<seqan3::sequence_file_input<>>{}(arguments.query_file);
    }

    arguments.treshold_was_set = parser.is_option_set("threshold");

    bool partitioned{false};
    seqan3::input_file_validator validator{};

    try
    {
        validator(arguments.index_file.string() + std::string{"_0"});
        partitioned = true;
    }
    catch (seqan3::validation_error const & e)
    {
        validator(arguments.index_file);
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
        std::ifstream is{partitioned ? arguments.index_file.string() + std::string{"_0"} :
                                       arguments.index_file.string(),
                         std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        raptor_index<> tmp{};
        tmp.load_parameters(iarchive);
        arguments.shape = tmp.shape();
        arguments.shape_size = arguments.shape.size();
        arguments.shape_weight = arguments.shape.count();
        arguments.window_size = tmp.window_size();
        arguments.parts = tmp.parts();
        arguments.compressed = tmp.compressed();
        arguments.bin_path = tmp.bin_path();
        if (arguments.is_socks)
            arguments.pattern_size = arguments.shape_size;
    }

    // ==========================================
    // Partitioned index: Check that all parts are available.
    // ==========================================
    if (partitioned)
    {
        for (size_t part{0}; part < arguments.parts; ++part)
        {
            validator(arguments.index_file.string() + std::string{"_"} + std::to_string(part));
        }
    }

    // ==========================================
    // Dispatch
    // ==========================================
    raptor_search(arguments);
}

} // namespace raptor

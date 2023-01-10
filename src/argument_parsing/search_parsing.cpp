// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <seqan3/io/views/async_input_buffer.hpp>

#include <raptor/argument_parsing/init_shared_meta.hpp>
#include <raptor/argument_parsing/search_parsing.hpp>
#include <raptor/argument_parsing/validators.hpp>
#include <raptor/dna4_traits.hpp>
#include <raptor/index.hpp>
#include <raptor/search/search.hpp>

namespace raptor
{

void init_search_parser(sharg::parser & parser, search_arguments & arguments)
{
    init_shared_meta(parser);
    parser.info.description.emplace_back("Searches a Raptor index using one or more sequences queries.");
    parser.info.examples = {
        "raptor search --error 2 --index raptor.index --query queries.fastq --output search.output"};

    parser.add_subsection("General options");
    parser.add_option(arguments.index_file,
                      sharg::config{.short_id = '\0',
                                    .long_id = "index",
                                    .description = arguments.is_socks
                                                     ? "Provide a valid path to an index."
                                                     : "Provide a valid path to an index. Parts: Without suffix _0",
                                    .required = true});
    parser.add_option(arguments.fpr,
                      sharg::config{.short_id = '\0',
                                    .long_id = "fpr",
                                    .description = "The false positive rate used for building the index.",
                                    .hidden = arguments.is_socks,
                                    .validator = sharg::arithmetic_range_validator{0, 1}});
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
    parser.add_option(arguments.threads,
                      sharg::config{.short_id = '\0',
                                    .long_id = "threads",
                                    .description = "The number of threads to use.",
                                    .validator = positive_integer_validator{}});

    parser.add_subsection("Threshold method options");
    parser.add_option(arguments.errors,
                      sharg::config{.short_id = '\0',
                                    .long_id = "error",
                                    .description = "The number of errors. Mutually exclusive with --threshold.",
                                    .hidden = arguments.is_socks,
                                    .validator = positive_integer_validator{true}});
    parser.add_option(arguments.threshold,
                      sharg::config{.short_id = '\0',
                                    .long_id = "threshold",
                                    .description = "Ratio of k-mers that need to be found for a hit to occur. "
                                                   "Mutually exclusive with --error.",
                                    .default_message = "None",
                                    .hidden = arguments.is_socks,
                                    .validator = sharg::arithmetic_range_validator{0, 1}});
    parser.add_option(arguments.pattern_size,
                      sharg::config{.short_id = '\0',
                                    .long_id = "query_length",
                                    .description = "The query length. Only used with --error.",
                                    .default_message = "Median of sequence lengths in query file",
                                    .hidden = arguments.is_socks});

    parser.add_subsection("Dynamic thresholding options");
    parser.add_line("\\fBThese option only take effect when using --error.\\fP");
    parser.add_option(arguments.tau,
                      sharg::config{.short_id = '\0',
                                    .long_id = "tau",
                                    .description = "The higher tau, the lower the threshold.",
                                    .hidden = arguments.is_socks,
                                    .validator = sharg::arithmetic_range_validator{0, 1}});
    parser.add_option(arguments.p_max,
                      sharg::config{.short_id = '\0',
                                    .long_id = "p_max",
                                    .description = "The higher p_max, the lower the threshold.",
                                    .hidden = arguments.is_socks,
                                    .validator = sharg::arithmetic_range_validator{0, 1}});
    parser.add_flag(
        arguments.cache_thresholds,
        sharg::config{
            .short_id = '\0',
            .long_id = "cache-thresholds",
            .description =
                "Stores the computed thresholds with an unique name next to the index. In the next search call "
                "using this option, the stored thresholds are re-used.\n"
                "Two files are stored:\n"
                "\\fBthreshold_*.bin\\fP: Depends on query_length, window, kmer/shape, errors, and tau.\n"
                "\\fBcorrection_*.bin\\fP: Depends on query_length, window, kmer/shape, p_max, and fpr."});

    parser.add_subsection("Advanced options", true);
    parser.add_flag(
        arguments.is_hibf,
        sharg::config{.short_id = '\0', .long_id = "hibf", .description = "Index is an HIBF.", .advanced = true});
    parser.add_flag(
        arguments.write_time,
        sharg::config{.short_id = '\0', .long_id = "time", .description = "Write timing file.", .advanced = true});
}

void search_parsing(sharg::parser & parser, bool const is_socks)
{
    search_arguments arguments{};
    arguments.is_socks = is_socks;
    init_search_parser(parser, arguments);
    parser.parse();

    // ==========================================
    // Various checks.
    // ==========================================

    if (parser.is_option_set("error") && parser.is_option_set("threshold"))
        throw sharg::parser_error{"You cannot set both error and threshold arguments."};

    std::filesystem::path output_directory = arguments.out_file.parent_path();
    std::error_code ec{};
    std::filesystem::create_directories(output_directory, ec);

    // GCOVR_EXCL_START
    if (!output_directory.empty() && ec)
        throw sharg::parser_error{
            seqan3::detail::to_string("Failed to create directory\"", output_directory.c_str(), "\": ", ec.message())};
    // GCOVR_EXCL_STOP

    if (!arguments.is_socks)
    {
        sharg::input_file_validator{raptor::detail::combined_extensions}(arguments.query_file);
    }

    if (std::filesystem::is_empty(arguments.query_file))
        throw sharg::parser_error{"The query file is empty."};

    bool partitioned{false};
    sharg::input_file_validator validator{};

    try
    {
        validator(arguments.index_file.string() + std::string{"_0"});
        partitioned = true;
    }
    catch (sharg::validation_error const & e)
    {
        validator(arguments.index_file);
    }

    // ==========================================
    // Process --query_length.
    // ==========================================
    if (!arguments.is_socks)
    {
        if (!parser.is_option_set("query_length"))
        {
            std::vector<uint64_t> sequence_lengths{};
            seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>> query_in{arguments.query_file};
            for (auto & [seq] : query_in | seqan3::views::async_input_buffer(16))
            {
                sequence_lengths.push_back(std::ranges::size(seq));
            }
            std::sort(sequence_lengths.begin(), sequence_lengths.end());
            arguments.pattern_size = sequence_lengths[sequence_lengths.size() / 2];
        }
    }

    // ==========================================
    // Read window and kmer size, and the bin paths.
    // ==========================================
    {
        std::ifstream is{partitioned ? arguments.index_file.string() + std::string{"_0"}
                                     : arguments.index_file.string(),
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

    if (arguments.pattern_size < arguments.window_size)
        throw sharg::parser_error{std::string{"The query size ("} + std::to_string(arguments.pattern_size)
                                  + ") is too short to be used with window size "
                                  + std::to_string(arguments.window_size) + '.'};

    // ==========================================
    // Temporary.
    // ==========================================
    if (arguments.shape_size != arguments.window_size && !parser.is_option_set("threshold")
        && !parser.is_option_set("fpr"))
    {
        std::cerr << "[WARNING] The search needs the FPR that was used for building the index.\n"
                  << "          Currently, the default value of " << std::setprecision(4) << arguments.fpr
                  << " is used.\n"
                  << "          If the index was built with a different FPR, the search results are not reliable.\n"
                  << "          The final version will store the FPR in the index and this parameter will be removed.\n"
                  << "          To disable this warning, explicitly pass the FPR to raptor search (--fpr 0.05).\n";
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

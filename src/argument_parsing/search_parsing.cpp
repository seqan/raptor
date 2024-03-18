// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::search_parsing.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

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
    parser.info.description.emplace_back("Queries a Raptor index.");
    parser.info.examples.emplace_back(
        "raptor search --index raptor.index --query queries.fastq --output search.output");
    parser.info.examples.emplace_back(
        "raptor search --index raptor.index --query queries.fastq --output search.output --threshold 0.7");
    parser.info.examples.emplace_back(
        "raptor search --index raptor.index --query queries.fastq --output search.output --error 2");
    parser.info.examples.emplace_back(
        "raptor search --index raptor.index --query queries.fastq --output search.output --error 2 --query_length 250");
    parser.info.synopsis.emplace_back("raptor search --index <file> --query <file> --output <file> [--threads "
                                      "<number>] [--quiet] [--error <number>|--threshold <number>] [--query_length "
                                      "<number>] [--tau <number>] [--pmax <number>] [--cache-thresholds]");
    parser.add_subsection("General options");
    parser.add_option(arguments.index_file,
                      sharg::config{.short_id = '\0',
                                    .long_id = "index",
                                    .description = "Provide a valid path to an index. Parts: Without suffix _0",
                                    .required = true});
    parser.add_option(arguments.query_file,
                      sharg::config{.short_id = '\0',
                                    .long_id = "query",
                                    .description = "Provide a path to the query file.",
                                    .required = true,
                                    .validator = sequence_file_validator{raptor::detail::combined_extensions}});
    parser.add_option(arguments.out_file,
                      sharg::config{.short_id = '\0',
                                    .long_id = "output",
                                    .description = "",
                                    .required = true,
                                    .validator = output_file_validator{}});
    parser.add_option(arguments.threads,
                      sharg::config{.short_id = '\0',
                                    .long_id = "threads",
                                    .description = "The number of threads to use.",
                                    .validator = positive_integer_validator{}});
    parser.add_flag(arguments.quiet,
                    sharg::config{.short_id = '\0',
                                  .long_id = "quiet",
                                  .description = "Do not print time and memory usage to stderr."});
    parser.add_option(arguments.timing_out,
                      sharg::config{.short_id = '\0',
                                    .long_id = "timing-output",
                                    .description = "Write time and memory usage to specified file (TSV format).",
                                    .validator = output_file_validator{}});

    parser.add_subsection("Threshold method options");
    parser.add_line("\\fBIf no option is set, --error " + std::to_string(arguments.errors)
                    + " will be used as default.\\fP");
    parser.add_line("Setting --query_length also skips validations associated with the query length, e.g., too long/"
                    "short queries or a high variance in the query lengths. Notably, queries that are shorter than the "
                    "window size result in undefined behavior.");
    parser.add_option(arguments.errors,
                      sharg::config{.short_id = '\0',
                                    .long_id = "error",
                                    .description = "The number of errors. Mutually exclusive with --threshold.",
                                    .validator = positive_integer_validator{true}});
    parser.add_option(arguments.threshold,
                      sharg::config{.short_id = '\0',
                                    .long_id = "threshold",
                                    .description = "Ratio of k-mers that need to be found for a hit to occur. "
                                                   "Mutually exclusive with --error.",
                                    .default_message = "None",
                                    .validator = sharg::arithmetic_range_validator{0, 1}});
    parser.add_option(arguments.query_length,
                      sharg::config{.short_id = '\0',
                                    .long_id = "query_length",
                                    .description =
                                        "The query length. Only influences the threshold when using --error. Enables "
                                        "skipping of the query length computation for both --error and --threshold.",
                                    .default_message = "Median of sequence lengths in query file"});

    parser.add_subsection("Dynamic thresholding options");
    parser.add_line("\\fBThese option have no effect when using --threshold or k-mer size == window size.\\fP");
    parser.add_line("A higher threshold means that less, in the extreme case even false negative, results will be "
                    "produced.");
    parser.add_line("A lower threshold means that more, possibly false positive, results will be produced.");
    parser.add_option(arguments.tau,
                      sharg::config{.short_id = '\0',
                                    .long_id = "tau",
                                    .description = "The higher tau, the lower the threshold.",
                                    .validator = sharg::arithmetic_range_validator{0, 1}});
    parser.add_option(arguments.p_max,
                      sharg::config{.short_id = '\0',
                                    .long_id = "p_max",
                                    .description = "The higher p_max, the lower the threshold.",
                                    .validator = sharg::arithmetic_range_validator{0, 1}});
    parser.add_flag(
        arguments.cache_thresholds,
        sharg::config{
            .short_id = '\0',
            .long_id = "cache-thresholds",
            .description =
                "Stores the computed thresholds with an unique name next to the index. In the next search call "
                "using this option, the stored thresholds are re-used. Two files are stored:"});
    parser.add_list_item("", "\\fBthreshold_*.bin\\fP: Depends on query_length, window, kmer/shape, errors, and tau.");
    parser.add_list_item("", "\\fBcorrection_*.bin\\fP: Depends on query_length, window, kmer/shape, p_max, and fpr.");
}

void search_parsing(sharg::parser & parser)
{
    search_arguments arguments{};
    arguments.wall_clock_timer.start();

    init_search_parser(parser, arguments);
    parser.parse();

    // ==========================================
    // Various checks.
    // ==========================================

    if (parser.is_option_set("error") && parser.is_option_set("threshold"))
        throw sharg::parser_error{"You cannot set both error and threshold arguments."};

    if (std::filesystem::is_empty(arguments.query_file))
        throw sharg::parser_error{"The query file is empty."};

    std::filesystem::path const partitioned_index_file = arguments.index_file.string() + "_0";
    bool const index_is_monolithic = std::filesystem::exists(arguments.index_file);
    bool const index_is_partitioned = std::filesystem::exists(partitioned_index_file);
    sharg::input_file_validator const index_validator{};

    if (index_is_monolithic && index_is_partitioned)
    {
        throw sharg::validation_error{sharg::detail::to_string("Ambiguous index. Both monolithic (",
                                                               arguments.index_file.c_str(),
                                                               ") and partitioned index (",
                                                               partitioned_index_file.c_str(),
                                                               ") exist. Please rename the monolithic index.")};
    }
    else if (index_is_partitioned)
    {
        index_validator(partitioned_index_file);
    }
    else
    {
        index_validator(arguments.index_file);
    }

    // ==========================================
    // Process --query_length.
    // ==========================================
    size_t min_query_length{arguments.query_length};
    size_t max_query_length{arguments.query_length};

    if (!parser.is_option_set("query_length"))
    {
        arguments.query_length_timer.start();
        std::vector<uint64_t> sequence_lengths{};
        seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>> query_in{arguments.query_file};

        for (auto && record : query_in | seqan3::views::async_input_buffer(1024))
            sequence_lengths.push_back(std::ranges::size(record.sequence()));

        std::ranges::sort(sequence_lengths);
        arguments.query_length = sequence_lengths[sequence_lengths.size() / 2];
        min_query_length = sequence_lengths.front();
        max_query_length = sequence_lengths.back();

        if (!parser.is_option_set("threshold") && max_query_length - min_query_length > arguments.query_length / 20u)
        {
            std::cerr << "[WARNING] There is variance in the provided queries. The shortest length is "
                      << min_query_length << ". The longest length is " << max_query_length
                      << ". The tresholding will use a single query length (" << arguments.query_length
                      << "). Therefore, results may be inprecise.\n";
        }
        arguments.query_length_timer.stop();
    }

    // We currently use counting_agent<uint16_t> and membership_agent (which uses uint16_t fixed).
    if (max_query_length > std::numeric_limits<uint16_t>::max())
    {
        std::cerr << "[WARNING] There are queries which exceed the maximum safely supported length of "
                  << std::numeric_limits<uint16_t>::max()
                  << ". Results may be wrong, especially when using window size == k-mer size. If you need longer "
                     "queries to be supported, please open an issue at https://github.com/seqan/raptor/issues.\n";
    }

    // ==========================================
    // Read window and kmer size, and the bin paths.
    // ==========================================
    {
        std::ifstream is{index_is_partitioned ? partitioned_index_file : arguments.index_file, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        raptor_index<> tmp{};
        tmp.load_parameters(iarchive);
        arguments.shape = tmp.shape();
        arguments.shape_size = arguments.shape.size();
        arguments.shape_weight = arguments.shape.count();
        arguments.window_size = tmp.window_size();
        arguments.parts = tmp.parts();
        arguments.bin_path = tmp.bin_path();
        arguments.fpr = tmp.fpr();
        arguments.is_hibf = tmp.is_hibf();
    }

    if (min_query_length < arguments.window_size)
        throw sharg::parser_error{sharg::detail::to_string("The (minimal) query length (",
                                                           min_query_length,
                                                           ") is too short to be used with window size ",
                                                           arguments.window_size,
                                                           '.')};

    // ==========================================
    // Partitioned index: Check that all parts are available.
    // ==========================================
    if (index_is_partitioned)
    {
        // GCOVR_EXCL_START
        std::string const index_path_base{[&partitioned_index_file]()
                                          {
                                              std::string_view sv = partitioned_index_file.c_str();
                                              assert(sv.size() > 0u);
                                              sv.remove_suffix(1u);
                                              return sv;
                                          }()};
        // GCOVR_EXCL_STOP
        for (size_t part{1u}; part < arguments.parts; ++part)
            index_validator(index_path_base + std::to_string(part));
    }

    // ==========================================
    // Dispatch
    // ==========================================
    raptor_search(arguments);

    arguments.wall_clock_timer.stop();
    if (!arguments.quiet)
        arguments.print_timings();
    if (parser.is_option_set("timing-output"))
        arguments.write_timings_to_file();
}

} // namespace raptor

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
                      "Provide a valid path to an IBF. Parts: Without suffix _0",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{});
    parser.add_option(arguments.query_file,
                      '\0',
                      "query",
                      "Provide a path to the query file.",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator<seqan3::sequence_file_input<>>{});
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
                      seqan3::option_spec::standard,
                      positive_integer_validator{true});
    parser.add_option(arguments.tau,
                      '\0',
                      "tau",
                      "Threshold for probabilistic models.",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{0, 1});
    parser.add_option(arguments.threshold,
                      '\0',
                      "threshold",
                      "If set, this threshold is used instead of the probabilistic models.",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{0, 1});
    parser.add_option(arguments.pattern_size,
                      '\0',
                      "pattern",
                      "Choose the pattern size. Default: Use median of sequence lengths in query file.");
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

void run_search(seqan3::argument_parser & parser)
{
    search_arguments arguments{};
    init_search_parser(parser, arguments);
    try_parsing(parser);

    // ==========================================
    // Various checks.
    // ==========================================

    if (parser.is_option_set("window"))
    {
        if (arguments.kmer_size > arguments.window_size)
            throw seqan3::argument_parser_error{"The k-mer size cannot be bigger than the window size."};
    }
    else
        arguments.window_size = arguments.kmer_size;

    arguments.treshold_was_set = parser.is_option_set("threshold");

    // ==========================================
    // Process --pattern.
    // ==========================================
    if (!arguments.pattern_size)
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
    // Dispatch
    // ==========================================
    raptor_search(arguments);
};

} // namespace raptor

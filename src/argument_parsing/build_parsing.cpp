// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Implements raptor::build_parsing.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <cereal/archives/json.hpp>

#include <chopper/configuration.hpp>

#include <raptor/argument_parsing/build_parsing.hpp>
#include <raptor/argument_parsing/compute_bin_size.hpp>
#include <raptor/argument_parsing/init_shared_meta.hpp>
#include <raptor/argument_parsing/parse_bin_path.hpp>
#include <raptor/argument_parsing/shared.hpp>
#include <raptor/argument_parsing/validators.hpp>
#include <raptor/build/hibf/bin_prefixes.hpp>
#include <raptor/build/raptor_build.hpp>

namespace raptor
{

inline bool read_chopper_config(chopper::configuration & config, std::filesystem::path const & config_file)
{
    std::ifstream file_in{config_file};
    std::string line{};
    std::stringstream config_stream{};

    std::getline(file_in, line, '\n');
    if (line != "##CONFIG:")
        return false;

    while (std::getline(file_in, line, '\n') && line.starts_with("##"))
    {
        if (line == "##ENDCONFIG")
            break;

        std::string_view sv{line};
        assert(line.size() > 1u);
        sv.remove_prefix(2u);
        config_stream << sv << '\n';
    }

    if (line != "##ENDCONFIG")
        return false; // GCOVR_EXCL_LINE

    cereal::JSONInputArchive iarchive(config_stream);
    iarchive(config);
    return true;
}

inline void parse_chopper_config(sharg::parser & parser, build_arguments & arguments)
{
    chopper::configuration config{};
    if (!read_chopper_config(config, arguments.bin_file))
        return;

    if (parser.is_option_set("kmer") && config.k != arguments.kmer_size)
    {
        std::cerr << sharg::detail::to_string(
            "[WARNING] Given k-mer size(",
            arguments.kmer_size,
            ") differs from k-mer size in the layout file (",
            config.k,
            "). The results may be suboptimal. If this was a conscious decision, you can ignore this warning.\n");
    }
    else
    {
        arguments.kmer_size = config.k;
    }

    validate_shape(parser, arguments);

    if (parser.is_option_set("hash") && config.num_hash_functions != arguments.hash)
    {
        std::cerr << sharg::detail::to_string(
            "[WARNING] Given hash function count (",
            arguments.hash,
            ") differs from hash function count in the layout file (",
            config.num_hash_functions,
            "). The results may be suboptimal. If this was a conscious decision, you can ignore this warning.\n");
    }
    else
        arguments.hash = config.num_hash_functions;

    if (parser.is_option_set("fpr") && config.false_positive_rate != arguments.fpr)
    {
        std::cerr << sharg::detail::to_string(
            "[WARNING] Given false positive rate (",
            arguments.fpr,
            ") differs from false positive rate in the layout file (",
            config.false_positive_rate,
            "). The results may be suboptimal. If this was a conscious decision, you can ignore this warning.");
    }
    else
        arguments.fpr = config.false_positive_rate;
}

inline void parse_shape_from_minimiser(sharg::parser & parser, build_arguments & arguments)
{
    if (parser.is_option_set("shape"))
        throw sharg::parser_error{"You cannot set --shape when using minimiser files as input."};
    if (parser.is_option_set("kmer"))
        throw sharg::parser_error{"You cannot set --kmer when using minimiser files as input."};
    if (parser.is_option_set("window"))
        throw sharg::parser_error{"You cannot set --window when using minimiser files as input."};

    std::filesystem::path header_file_path = arguments.bin_path[0][0];
    header_file_path.replace_extension("header");
    std::ifstream file_stream{header_file_path};
    std::string shape_string{};
    file_stream >> shape_string >> arguments.window_size;

    uint64_t tmp{};
    std::from_chars(shape_string.data(), shape_string.data() + shape_string.size(), tmp, 2);

    arguments.shape = seqan3::shape{seqan3::bin_literal{tmp}};
}

void init_build_parser(sharg::parser & parser, build_arguments & arguments)
{
    init_shared_meta(parser);
    parser.info.description.emplace_back("Constructs a Raptor index.");
    parser.info.description.emplace_back("The input may be a layout file from \\fBraptor layout\\fP, a list of "
                                         "minimiser files produced from \\fBraptor prepare\\fP, or a file with a list "
                                         "of files to process.");
    parser.info.examples.emplace_back("raptor build --kmer 19 --window 23 --fpr 0.05 --output raptor.index bins.list");
    parser.info.examples.emplace_back("raptor build --shape 11011 --window 8 --output raptor.index bins.list");
    parser.info.examples.emplace_back(
        "raptor build --kmer 32 --window 32 --hash 3 --parts 4 --output raptor.index bins.list");
    parser.info.examples.emplace_back("raptor build --fpr 0.05 --output raptor.index minimiser.list");
    parser.info.examples.emplace_back("raptor build --output raptor.index raptor.layout");
    parser.info.examples.emplace_back("raptor build --fpr 0.05 --output raptor.index raptor.layout");
    parser.info.synopsis.emplace_back("raptor build --output <file> [--threads <number>] [--verbose] [--kmer <number>"
                                      "|--shape <01-pattern>] [--window <number>] [--fpr <number>] [--hash <number>] "
                                      "[--parts <number>] [--compressed] [--] <INPUT>");

    parser.add_positional_option(
        arguments.bin_file,
        sharg::config{.description = "A layout file from \\fBraptor layout\\fP, or a file containing file names. "
                                   + bin_validator{}.get_help_page_message(),
                      .validator = sharg::input_file_validator{}});

    parser.add_subsection("General options");
    parser.add_option(arguments.out_path,
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
    parser.add_flag(
        arguments.verbose,
        sharg::config{.short_id = '\0', .long_id = "verbose", .description = "Print time and memory usage."});

    parser.add_subsection("k-mer options");
    parser.add_option(
        arguments.kmer_size,
        sharg::config{.short_id = '\0',
                      .long_id = "kmer",
                      .description = "The k-mer size.",
                      .default_message = std::to_string(arguments.kmer_size) + ", or read from layout file",
                      .validator = sharg::arithmetic_range_validator{1, 32}});
    parser.add_option(arguments.window_size,
                      sharg::config{.short_id = '\0',
                                    .long_id = "window",
                                    .description = "The window size.",
                                    .default_message = "k-mer size",
                                    .validator = positive_integer_validator{}});
    parser.add_option(
        arguments.shape_string,
        sharg::config{.short_id = '\0',
                      .long_id = "shape",
                      .description =
                          "The shape to use for k-mers. Mutually exclusive with --kmer. Parsed from right to left.",
                      .default_message = "11111111111111111111 (a k-mer of size 20), or read from layout file",
                      .validator = sharg::regex_validator{"[01]+"}});

    parser.add_subsection("Index options");
    parser.add_option(arguments.fpr,
                      sharg::config{.short_id = '\0',
                                    .long_id = "fpr",
                                    .description = "The false positive rate.",
                                    .default_message = std::to_string(arguments.fpr) + ", or read from layout file",
                                    .validator = sharg::arithmetic_range_validator{0.0, 1.0}});
    parser.add_option(arguments.hash,
                      sharg::config{.short_id = '\0',
                                    .long_id = "hash",
                                    .description = "The number of hash functions to use.",
                                    .default_message = std::to_string(arguments.hash) + ", or read from layout file",
                                    .validator = sharg::arithmetic_range_validator{1, 5}});
    parser.add_option(arguments.parts,
                      sharg::config{.short_id = '\0',
                                    .long_id = "parts",
                                    .description = "Splits the index in this many parts. Not available for the HIBF.",
                                    .validator = power_of_two_validator{}});
    parser.add_flag(
        arguments.compressed,
        sharg::config{.short_id = '\0', .long_id = "compressed", .description = "Build a compressed index."});
}

bool input_is_pack_file(std::filesystem::path const & path)
{
    std::ifstream file{path};
    std::string line{};
    while (std::getline(file, line) && line.starts_with("##")) // Skip parameter information
    {}
    return line.starts_with(raptor::hibf::pack_file_first_line_prefix);
}

void build_parsing(sharg::parser & parser)
{
    build_arguments arguments{};
    arguments.wall_clock_timer.start();

    init_build_parser(parser, arguments);
    parser.parse();

    if (std::filesystem::is_empty(arguments.bin_file))
        throw sharg::parser_error{"The input file is empty."};

    arguments.is_hibf = input_is_pack_file(arguments.bin_file);

    if (arguments.is_hibf && arguments.parts != 1u)
        throw sharg::parser_error{"The HIBF cannot yet be partitioned."};

    parse_bin_path(arguments);

    if (arguments.is_hibf)
        parse_chopper_config(parser, arguments);

    if (!arguments.input_is_minimiser)
        validate_shape(parser, arguments);
    else
        parse_shape_from_minimiser(parser, arguments);

    if (!arguments.is_hibf && arguments.parts == 1u)
        arguments.bits = compute_bin_size(arguments);

    raptor_build(arguments);

    arguments.wall_clock_timer.stop();
    arguments.print_timings();
}

} // namespace raptor

// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <raptor/argument_parsing/build_parsing.hpp>
#include <raptor/argument_parsing/compute_bin_size.hpp>
#include <raptor/argument_parsing/init_shared_meta.hpp>
#include <raptor/argument_parsing/parse_bin_path.hpp>
#include <raptor/argument_parsing/shared.hpp>
#include <raptor/argument_parsing/validators.hpp>
#include <raptor/build/raptor_build.hpp>

namespace raptor
{

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
    parser.info.description.emplace_back("Constructs a Raptor index of large collections of nucleotide sequences.");
    parser.info.examples = {"raptor build --kmer 19 --window 23 --fpr 0.05 --output raptor.index all_bin_paths.txt"};

    parser.add_positional_option(
        arguments.bin_file,
        sharg::config{.description = "File containing file names. " + bin_validator{}.get_help_page_message(),
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

    parser.add_subsection("k-mer options");
    parser.add_option(arguments.kmer_size,
                      sharg::config{.short_id = '\0',
                                    .long_id = "kmer",
                                    .description = "The k-mer size.",
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
                      .default_message = "11111111111111111111 (a k-mer of size 20)",
                      .validator = sharg::regex_validator{"[01]+"}});

    parser.add_subsection("Index options");
    parser.add_option(arguments.fpr,
                      sharg::config{.short_id = '\0',
                                    .long_id = "fpr",
                                    .description = "The false positive rate.",
                                    .validator = sharg::arithmetic_range_validator{0.0, 1.0}});
    parser.add_option(arguments.hash,
                      sharg::config{.short_id = '\0',
                                    .long_id = "hash",
                                    .description = "The number of hash functions to use.",
                                    .validator = sharg::arithmetic_range_validator{1, 5}});
    parser.add_option(arguments.parts,
                      sharg::config{.short_id = '\0',
                                    .long_id = "parts",
                                    .description = "Splits the index in this many parts. Not available for the HIBF.",
                                    .validator = power_of_two_validator{}});
    parser.add_flag(
        arguments.compressed,
        sharg::config{.short_id = '\0', .long_id = "compressed", .description = "Build a compressed index."});
    parser.add_flag(arguments.is_hibf,
                    sharg::config{.short_id = '\0', .long_id = "hibf", .description = "Build an HIBF."});
}

void build_parsing(sharg::parser & parser)
{
    build_arguments arguments{};
    init_build_parser(parser, arguments);
    parser.parse();

    if (arguments.is_hibf && arguments.parts != 1u)
        throw sharg::parser_error{"The HIBF cannot yet be partitioned."};

    parse_bin_path(arguments);

    if (!arguments.input_is_minimiser)
        validate_shape(parser, arguments);
    else
        parse_shape_from_minimiser(parser, arguments);

    if (!arguments.is_hibf && arguments.parts == 1u)
        arguments.bits = compute_bin_size(arguments);

    raptor_build(arguments);
}

} // namespace raptor

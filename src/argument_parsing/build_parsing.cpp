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
#include <raptor/build/raptor_build.hpp>

#include <hibf/layout/prefixes.hpp>

namespace raptor
{

inline void read_chopper_config(chopper::configuration & config, std::filesystem::path const & config_file)
{
    std::ifstream file_in{config_file};
    config.read_from(file_in);
}

inline void parse_chopper_config(sharg::parser & parser, build_arguments & arguments)
{
    if (parser.is_option_set("hash"))
        throw sharg::parser_error{"You cannot set --hash when using a layout file as input."};
    if (parser.is_option_set("fpr"))
        throw sharg::parser_error{"You cannot set --fpr when using a layout file as input."};

    chopper::configuration config{};

    read_chopper_config(config, arguments.bin_file);

    if (parser.is_option_set("kmer") && config.k != arguments.kmer_size)
    {
        std::cerr << sharg::detail::to_string(
            "[WARNING] Given k-mer size (",
            arguments.kmer_size,
            ") differs from k-mer size in the layout file (",
            config.k,
            "). The results may be suboptimal. If this was a conscious decision, you can ignore this warning.\n");
    }
    else if (parser.is_option_set("shape") && config.k != arguments.shape_string.size())
    {
        std::cerr << sharg::detail::to_string(
            "[WARNING] Given size of the shape (",
            arguments.shape_string.size(),
            ") differs from k-mer size in the layout file (",
            config.k,
            "). The results may be suboptimal. If this was a conscious decision, you can ignore this warning.\n");
    }
    else
    {
        arguments.kmer_size = config.k;
    }

    validate_shape(parser, arguments);
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
    parser.info.examples.emplace_back("raptor build --input bins.list --kmer 19 --window 23 --fpr 0.05 --output "
                                      "raptor.index");
    parser.info.examples.emplace_back("raptor build --input bins.list --shape 11011 --window 8 --output raptor.index");
    parser.info.examples.emplace_back("raptor build --input bins.list --kmer 32 --window 32 --hash 3 --parts 4 "
                                      "--output raptor.index");
    parser.info.examples.emplace_back("raptor build --input minimiser.list --fpr 0.05 --output raptor.index");
    parser.info.examples.emplace_back("raptor build --input raptor.layout --output raptor.index");
    parser.info.synopsis.emplace_back("raptor build --input <file> --output <file> [--threads <number>] [--quiet] "
                                      "[--kmer <number>|--shape <01-pattern>] [--window <number>] [--fpr <number>] "
                                      "[--hash <number>] [--parts <number>]");

    parser.add_subsection("General options");
    parser.add_option(
        arguments.bin_file,
        sharg::config{.short_id = '\0',
                      .long_id = "input",
                      .description = "A layout file from \\fBraptor layout\\fP, or a file containing file names. "
                                   + bin_validator{}.get_help_page_message(),
                      .required = true,
                      .validator = sharg::input_file_validator{}});
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
    parser.add_flag(arguments.quiet,
                    sharg::config{.short_id = '\0',
                                  .long_id = "quiet",
                                  .description = "Do not print time and memory usage to stderr."});
    parser.add_option(arguments.timing_out,
                      sharg::config{.short_id = '\0',
                                    .long_id = "timing-output",
                                    .description = "Write time and memory usage to specified file (TSV format).",
                                    .validator = output_file_validator{}});

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

    // GCOVR_EXCL_START
    // Adding additional cwl information that currently aren't supported by sharg and tdl.
    tdl::post_process_cwl = [](YAML::Node & node)
    {
        node["requirements"] = YAML::Load(R"-(
                                              InlineJavascriptRequirement: {}
                                              InitialWorkDirRequirement:
                                                listing:
                                                  - $(inputs.sequences.flat())
                                                  - entryname: input_bins_filepaths.txt
                                                    entry: |
                                                      ${
                                                         var bins = "";
                                                         for (var i = 0; i < inputs.sequences.length; i++) {
                                                            var currentBin = inputs.sequences[i];
                                                            for (var j = 0; j < currentBin.length; j++) {
                                                              bins += currentBin[j].basename + " ";
                                                            }
                                                            bins += "\n";
                                                         }
                                                         return bins;
                                                      }
                                            )-");
        auto inputs = node["inputs"];
        inputs.remove("input");
        inputs["output_name"] = inputs["output"];
        inputs.remove("output");
        inputs["sequences"] = YAML::Load(R"-(
                                             type:
                                               type: array
                                               items:
                                                 type: array
                                                 items: File
                                            )-");
        node["outputs"] = YAML::Load(R"-(
                                         index:
                                           type: File
                                           outputBinding:
                                             glob: $(inputs.output_name)
                                       )-");
        node["arguments"] = YAML::Load(R"-(
                                            - --input
                                            - input_bins_filepaths.txt
                                          )-");
    };
    // GCOVR_EXCL_STOP
}

bool input_is_pack_file(std::filesystem::path const & path)
{
    char const first_char = std::ifstream{path}.get();

    if (first_char == '#')
        throw sharg::parser_error{"The input file was determined to be a layout. However, the first line starts with "
                                  "'#' (old format) instead of '@' (new format)."};

    return seqan::hibf::prefix::meta_header.front() == first_char;
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
    if (!arguments.quiet)
        arguments.print_timings();
    if (parser.is_option_set("timing-output"))
        arguments.write_timings_to_file();
}

} // namespace raptor

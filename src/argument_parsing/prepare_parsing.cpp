// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <raptor/argument_parsing/init_shared_meta.hpp>
#include <raptor/argument_parsing/parse_bin_path.hpp>
#include <raptor/argument_parsing/prepare_arguments.hpp>
#include <raptor/argument_parsing/prepare_parsing.hpp>
#include <raptor/argument_parsing/shared.hpp>
#include <raptor/argument_parsing/validators.hpp>
#include <raptor/prepare/compute_minimiser.hpp>

namespace raptor
{

void init_prepare_parser(sharg::parser & parser, prepare_arguments & arguments)
{
    init_shared_meta(parser);
    parser.info.description.emplace_back("Computes minimisers for the use with raptor build. Creates minimiser and "
                                         "header files for each given file in the input file.");
    // parser.info.examples = {};

    parser.add_positional_option(
        arguments.bin_file,
        sharg::config{.description = "File containing file names. " + bin_validator{}.get_help_page_message(),
                      .validator = sharg::input_file_validator{}});

    parser.add_subsection("General options");
    parser.add_option(arguments.out_dir,
                      sharg::config{.short_id = '\0',
                                    .long_id = "output",
                                    .description = "",
                                    .required = true,
                                    .validator = output_directory_validator{}});
    parser.add_option(arguments.threads,
                      sharg::config{.short_id = '\0',
                                    .long_id = "threads",
                                    .description = "The number of threads to use.",
                                    .validator = positive_integer_validator{}});
    parser.add_flag(
        arguments.verbose,
        sharg::config{.short_id = '\0', .long_id = "verbose", .description = "Print time and memory usage."});

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

    parser.add_subsection("Processing options");
    parser.add_option(arguments.kmer_count_cutoff,
                      sharg::config{.short_id = '\0',
                                    .long_id = "kmer-count-cutoff",
                                    .description = "Only store k-mers with at least (>=) x occurrences. "
                                                   "Mutually exclusive with --use-filesize-dependent-cutoff.",
                                    .validator = sharg::arithmetic_range_validator{1, 254}});
    parser.add_flag(arguments.use_filesize_dependent_cutoff,
                    sharg::config{.short_id = '\0',
                                  .long_id = "use-filesize-dependent-cutoff",
                                  .description = "Apply cutoffs from Mantis(Pandey et al., 2018). "
                                                 "Mutually exclusive with --kmer-count-cutoff."});
}

void prepare_parsing(sharg::parser & parser)
{
    prepare_arguments arguments{};
    arguments.wall_clock_timer.start();

    init_prepare_parser(parser, arguments);
    parser.parse();

    if (parser.is_option_set("kmer-count-cutoff") && parser.is_option_set("use-filesize-dependent-cutoff"))
        throw sharg::parser_error{"You cannot use both --kmer-count-cutoff and --use-filesize-dependent-cutoff."};

    validate_shape(parser, arguments);

    parse_bin_path(arguments);

    compute_minimiser(arguments);

    arguments.wall_clock_timer.stop();
    arguments.print_timings();
}

} // namespace raptor

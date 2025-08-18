// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::prepare_parsing.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

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
    parser.info.short_description = "Computes minimisers for the use with raptor layout and raptor build";
    parser.info.description.emplace_back(
        "Computes minimisers for the use with \\fBraptor layout\\fP and \\fBraptor build\\fP.");
    parser.info.description.emplace_back("Can continue where it left off after a crash or in multiple runs.");
    parser.info.examples.emplace_back("raptor prepare --input bins.list --output some_directory --kmer 20 --window 24");
    parser.info.examples.emplace_back("raptor prepare --input bins.list --output some_directory --kmer-count-cutoff 2");
    parser.info.examples.emplace_back("raptor prepare --input bins.list --output some_directory "
                                      "--use-filesize-dependent-cutoff");
    parser.info.synopsis.emplace_back(
        "raptor prepare --input <file> --output <directory> [--threads <number>] [--quiet] [--kmer <number>|--shape "
        "<01-pattern>] [--window <number>] [--kmer-count-cutoff <number>|--use-filesize-dependent-cutoff]");

    parser.add_subsection("General options");
    parser.add_option(
        arguments.bin_file,
        sharg::config{.short_id = '\0',
                      .long_id = "input",
                      .description = "File containing file names. " + bin_validator{}.get_help_page_message(),
                      .required = true,
                      .validator = sharg::input_file_validator{}});
    parser.add_option(arguments.out_dir,
                      sharg::config{.short_id = '\0',
                                    .long_id = "output",
                                    .description = "",
                                    .required = true,
                                    .validator = output_directory_validator{}});
    parser.add_list_item("",
                         "Will create a \\fBminimiser.list\\fP inside the output directory. This file contains a "
                         "list of generated minimiser files, in the same order as the input.");
    parser.add_list_item("",
                         "\\fBWhen you manually delete a .in_progress file, also delete the corresponding .header and "
                         ".minimiser file!\\fP");
    parser.add_list_item("", "Created output files for each file:");
    parser.add_list_item("", "\\fB*.header\\fP: Contains the shape, window size, cutoff and minimiser count.");
    parser.add_list_item("", "\\fB*.minimiser\\fP: Contains binary minimiser values, one minimiser per line.");
    parser.add_list_item(
        "",
        "\\fB*.in_progress\\fP: Temporary file to track process. Deleted after finishing computation.");
    parser.add_option(arguments.threads,
                      sharg::config{.short_id = '\0',
                                    .long_id = "threads",
                                    .description = "The number of threads to use.",
                                    .validator = positive_integer_validator{}});
    parser.add_flag(
        arguments.quiet,
        sharg::config{.short_id = '\0', .long_id = "quiet", .description = "Do not print time and memory usage."});

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

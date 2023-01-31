// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <raptor/argument_parsing/init_shared_meta.hpp>
#include <raptor/argument_parsing/parse_bin_path.hpp>
#include <raptor/argument_parsing/upgrade_parsing.hpp>
#include <raptor/argument_parsing/validators.hpp>
#include <raptor/upgrade/upgrade.hpp>

namespace raptor
{

void init_upgrade_parser(sharg::parser & parser, upgrade_arguments & arguments)
{
    init_shared_meta(parser);
    parser.info.description.emplace_back("Upgrades a Raptor index created with an older version of Raptor to be"
                                         " compatible with the currently used version of Raptor.");

    parser.add_option(arguments.bin_file,
                      sharg::config{.short_id = '\0',
                                    .long_id = "bins",
                                    .description = "File containing one file per line per bin.",
                                    .required = true,
                                    .validator = sharg::input_file_validator{}});
    parser.add_option(arguments.in_file,
                      sharg::config{.short_id = '\0',
                                    .long_id = "input",
                                    .description = "The index to upgrade. Parts: Without suffix _0",
                                    .required = true});
    parser.add_option(
        arguments.out_file,
        sharg::config{.short_id = '\0', .long_id = "output", .description = "Path to new index.", .required = true});
    parser.add_option(arguments.window_size,
                      sharg::config{.short_id = '\0',
                                    .long_id = "window",
                                    .description = "The original window size.",
                                    .required = true,
                                    .validator = positive_integer_validator{}});
    parser.add_option(arguments.kmer_size,
                      sharg::config{.short_id = '\0',
                                    .long_id = "kmer",
                                    .description = "The original kmer size.",
                                    .required = true,
                                    .validator = sharg::arithmetic_range_validator{1, 32}});
    parser.add_option(arguments.parts,
                      sharg::config{.short_id = '\0',
                                    .long_id = "parts",
                                    .description = "Original index consisted of this many parts.",
                                    .validator = power_of_two_validator{}});

    parser.add_flag(
        arguments.compressed,
        sharg::config{.short_id = '\0', .long_id = "compressed", .description = "Original index was compressed."});
}

void upgrade_parsing(sharg::parser & parser)
{
    upgrade_arguments arguments{};
    init_upgrade_parser(parser, arguments);
    parser.parse();

    // ==========================================
    // Various checks.
    // ==========================================
    if (arguments.kmer_size > arguments.window_size)
        throw sharg::parser_error{"The k-mer size cannot be bigger than the window size."};

    arguments.shape = seqan3::shape{seqan3::ungapped{arguments.kmer_size}};

    std::filesystem::path output_directory = arguments.out_file.parent_path();
    std::error_code ec{};
    std::filesystem::create_directories(output_directory, ec);

    // GCOVR_EXCL_START
    if (!output_directory.empty() && ec)
        throw sharg::parser_error{
            seqan3::detail::to_string("Failed to create directory\"", output_directory.c_str(), "\": ", ec.message())};
    // GCOVR_EXCL_STOP

    if (arguments.parts == 1)
    {
        sharg::input_file_validator{}(arguments.in_file);
    }
    else
    {
        sharg::input_file_validator validator{};
        for (size_t part{0}; part < arguments.parts; ++part)
        {
            validator(arguments.in_file.string() + std::string{"_"} + std::to_string(part));
        }
    }

    // ==========================================
    // Process bin_path
    // ==========================================
    parse_bin_path(arguments);

    // ==========================================
    // Dispatch
    // ==========================================
    raptor_upgrade(arguments);
}

} // namespace raptor

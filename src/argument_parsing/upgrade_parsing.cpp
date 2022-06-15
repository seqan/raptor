// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <raptor/argument_parsing/init_shared_meta.hpp>
#include <raptor/argument_parsing/parse_bin_path.hpp>
#include <raptor/argument_parsing/upgrade_parsing.hpp>
#include <raptor/argument_parsing/validators.hpp>
#include <raptor/upgrade/upgrade.hpp>

namespace raptor
{

void init_upgrade_parser(seqan3::argument_parser & parser, upgrade_arguments & arguments)
{
    init_shared_meta(parser);
    parser.add_option(arguments.bin_file,
                      '\0',
                      "bins",
                      "File containing one file per line per bin.",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{});
    parser.add_option(arguments.in_file,
                      '\0',
                      "input",
                      "The index to upgrade. Parts: Without suffix _0",
                      seqan3::option_spec::required);
    // clang-format off
    parser.add_option(arguments.out_file,
                      '\0',
                      "output",
                      "Path to new index.",
                      seqan3::option_spec::required);
    // clang-format on
    parser.add_option(arguments.window_size,
                      '\0',
                      "window",
                      "The original window size.",
                      seqan3::option_spec::required,
                      positive_integer_validator{});
    parser.add_option(arguments.kmer_size,
                      '\0',
                      "kmer",
                      "The original kmer size.",
                      seqan3::option_spec::required,
                      seqan3::arithmetic_range_validator{1, 32});
    parser.add_option(arguments.parts,
                      '\0',
                      "parts",
                      "Original index consisted of this many parts.",
                      seqan3::option_spec::standard,
                      power_of_two_validator{});
    // clang-format off
    parser.add_flag(arguments.compressed,
                    '\0',
                    "compressed",
                    "Original index was compressed.");
    // clang-format on
}

void upgrade_parsing(seqan3::argument_parser & parser)
{
    upgrade_arguments arguments{};
    init_upgrade_parser(parser, arguments);
    parser.parse();

    // ==========================================
    // Various checks.
    // ==========================================
    if (arguments.kmer_size > arguments.window_size)
        throw seqan3::argument_parser_error{"The k-mer size cannot be bigger than the window size."};

    arguments.shape = seqan3::shape{seqan3::ungapped{arguments.kmer_size}};

    std::filesystem::path output_directory = arguments.out_file.parent_path();
    std::error_code ec{};
    std::filesystem::create_directories(output_directory, ec);

    // GCOVR_EXCL_START
    if (!output_directory.empty() && ec)
        throw seqan3::argument_parser_error{
            seqan3::detail::to_string("Failed to create directory\"", output_directory.c_str(), "\": ", ec.message())};
    // GCOVR_EXCL_STOP

    if (arguments.parts == 1)
    {
        seqan3::input_file_validator{}(arguments.in_file);
    }
    else
    {
        seqan3::input_file_validator validator{};
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

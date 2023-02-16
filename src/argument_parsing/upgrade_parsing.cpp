// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Implements raptor::upgrade_parsing.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <raptor/argument_parsing/init_shared_meta.hpp>
#include <raptor/argument_parsing/parse_bin_path.hpp>
#include <raptor/argument_parsing/upgrade_parsing.hpp>
#include <raptor/argument_parsing/validators.hpp>
#include <raptor/build/partition_config.hpp>
#include <raptor/index.hpp>
#include <raptor/upgrade/upgrade.hpp>

namespace raptor
{

void init_upgrade_parser(sharg::parser & parser, upgrade_arguments & arguments)
{
    init_shared_meta(parser);
    parser.info.description.emplace_back("Upgrades a Raptor index created with Raptor 2.0 to be"
                                         " compatible with Raptor 3.0.");
    parser.info.description.emplace_back("The only new parameter need is the false positive rate. The false positive"
                                         " rate affects the search results. This can be done in three different ways:");
    parser.info.description.emplace_back("\\fB1)\\fP Pass the false positive rate via --fpr.");
    parser.info.description.emplace_back("\\fB2)\\fP The false positive rate can be automatically determined if the"
                                         " paths of the files used to build the index are still available.");
    parser.info.description.emplace_back("\\fB3)\\fP Pass a file containing the path to the original files, one line"
                                         " per file. The false positive rate can then be automatically determined. The"
                                         " order of the files does not matter. The file with the most k-mers will"
                                         " determine the false positive rate.");
    parser.info.examples.emplace_back("raptor upgrade --input old.index --output new.index");
    parser.info.examples.emplace_back("raptor upgrade --input old.index --output new.index --fpr 0.05");
    parser.info.examples.emplace_back("raptor upgrade --input old.index --output new.index --bins bins.list");
    parser.info.synopsis.emplace_back("raptor upgrade --input <file> --output <file> [--fpr <number>|--bins <file>]");

    parser.add_option(arguments.fpr,
                      sharg::config{.short_id = '\0',
                                    .long_id = "fpr",
                                    .description = "The false positive rate. Mutually exclusive with --bins.",
                                    .default_message = "None",
                                    .validator = sharg::arithmetic_range_validator{0.0, 1.0}});
    parser.add_option(
        arguments.bin_file,
        sharg::config{.short_id = '\0',
                      .long_id = "bins",
                      .description = "File containing one file per line per bin. Mutually exclusive with --fpr.",
                      .default_message = "None",
                      .required = false,
                      .validator = sharg::input_file_validator{}});
    parser.add_option(arguments.index_file,
                      sharg::config{.short_id = '\0',
                                    .long_id = "input",
                                    .description = "The index to upgrade. Parts: Without suffix _0",
                                    .required = true});
    parser.add_option(
        arguments.output_file,
        sharg::config{.short_id = '\0', .long_id = "output", .description = "Path to new index.", .required = true});
}

void upgrade_parsing(sharg::parser & parser)
{
    upgrade_arguments arguments{};
    init_upgrade_parser(parser, arguments);
    parser.parse();

    if (parser.is_option_set("fpr") && parser.is_option_set("bins"))
        throw sharg::validation_error{"You cannot set both --fpr and --bins."};

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

    if (parser.is_option_set("bins"))
        parse_bin_path(arguments);

    {
        std::ifstream is{index_is_partitioned ? partitioned_index_file : arguments.index_file, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        raptor_index<> tmp{};
        tmp.load_old_parameters(iarchive);
        arguments.shape = tmp.shape();
        arguments.window_size = tmp.window_size();
        arguments.parts = tmp.parts();
        arguments.compressed = tmp.compressed();
        if (arguments.bin_path.empty() && !parser.is_option_set("fpr"))
        {
            arguments.bin_path = tmp.bin_path();
            bin_validator{}(arguments.bin_path);
            arguments.input_is_minimiser = arguments.bin_path[0][0].ends_with(".minimiser");
        }
    }

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

    raptor_upgrade(arguments);
}

} // namespace raptor

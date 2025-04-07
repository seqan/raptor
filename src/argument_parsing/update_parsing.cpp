// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::update_parsing.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <algorithm>   // for find_if
#include <charconv>    // for from_chars
#include <filesystem>  // for operator<<, operator>>, path
#include <fstream>     // for operator<<, basic_ifstream, basic_ios, basic_istream
#include <iomanip>     // for operator<<, quoted
#include <string>      // for operator+, basic_string, char_traits, to_string
#include <string_view> // for basic_string_view, operator==, string_view
#include <vector>      // for vector

#include <cereal/archives/binary.hpp> // for BinaryInputArchive
#include <cereal/specialize.hpp>      // for specialization

#include <sharg/auxiliary.hpp>        // for parser_meta_data
#include <sharg/config.hpp>           // for config
#include <sharg/detail/to_string.hpp> // for to_string
#include <sharg/parser.hpp>           // for parser
#include <sharg/validators.hpp>       // for input_file_validator, output_file_open_options

#include <seqan3/search/kmer_index/shape.hpp> // for shape

#include <hibf/interleaved_bloom_filter.hpp> // for interleaved_bloom_filter

#include <raptor/argument_parsing/init_shared_meta.hpp> // for init_shared_meta
#include <raptor/argument_parsing/update_arguments.hpp> // for update_arguments
#include <raptor/argument_parsing/update_parsing.hpp>   // for update_parsing
#include <raptor/argument_parsing/validators.hpp>       // for output_file_validator, positive_integer_validator
#include <raptor/index.hpp>                             // for raptor_index
#include <raptor/update/update.hpp>                     // for raptor_update

namespace raptor
{

void init_update_parser(sharg::parser & parser, update_arguments & arguments)
{
    init_shared_meta(parser);
    parser.info.description.emplace_back("Updates a Raptor index.");

    parser.add_subsection("General options");
    parser.add_option(arguments.index_file,
                      sharg::config{.short_id = '\0',
                                    .long_id = "index",
                                    .description = "Path to an index.",
                                    .required = true,
                                    .validator = sharg::input_file_validator{}});
    parser.add_option(arguments.out_path,
                      sharg::config{.short_id = '\0',
                                    .long_id = "output",
                                    .description = "",
                                    .required = true,
                                    .validator = output_file_validator{sharg::output_file_open_options::create_new}});
    parser.add_option(arguments.threads,
                      sharg::config{.short_id = '\0',
                                    .long_id = "threads",
                                    .description = "The number of threads to use.",
                                    .validator = positive_integer_validator{}});
}

void init_delete_parser(sharg::parser & parser, update_arguments & arguments)
{
    init_update_parser(parser, arguments);
    parser.add_subsection("Deletion options");
    parser.add_option(arguments.user_bins_to_delete,
                      sharg::config{
                          .short_id = '\0',
                          .long_id = "delete",
                          .description = "UB to delete",
                          .required = true,
                      });
}

void init_insert_parser(sharg::parser & parser, update_arguments & arguments)
{
    init_update_parser(parser, arguments);
    parser.add_subsection("Insertion options");
    parser.add_option(arguments.user_bin_to_insert,
                      sharg::config{.short_id = '\0',
                                    .long_id = "insert",
                                    .description = "UB to insert",
                                    .required = true,
                                    .validator = sharg::input_file_validator{}});
}

void update_parsing(sharg::parser & parser)
{
    init_shared_meta(parser);
    parser.info.description.emplace_back("Updates a Raptor index.");
    parser.add_subcommands({"delete", "insert"});
    parser.parse();

    sharg::parser & sub_parser = parser.get_sub_parser();
    update_arguments arguments{};

    if (sub_parser.info.app_name == std::string_view{"Raptor-update-delete"})
        init_delete_parser(sub_parser, arguments);
    if (sub_parser.info.app_name == std::string_view{"Raptor-update-insert"})
        init_insert_parser(sub_parser, arguments);

    sub_parser.parse();

    // ==========================================
    // Read window and kmer size, and the bin paths.
    // ==========================================
    {
        std::ifstream is{arguments.index_file, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        raptor_index<> tmp{};
        tmp.load_parameters(iarchive);
        arguments.shape = tmp.shape();
        arguments.shape_size = arguments.shape.size();
        arguments.shape_weight = arguments.shape.count();
        arguments.window_size = tmp.window_size();
        arguments.parts = tmp.parts();
        // arguments.bin_path = tmp.bin_path();
        arguments.fpr = tmp.fpr();
        arguments.is_hibf = tmp.is_hibf();
    }

    raptor_update(arguments);
}

} // namespace raptor

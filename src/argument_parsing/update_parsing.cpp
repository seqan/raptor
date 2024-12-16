// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::update_parsing.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <raptor/argument_parsing/init_shared_meta.hpp>
#include <raptor/argument_parsing/parse_bin_path.hpp>
#include <raptor/argument_parsing/update_parsing.hpp>
#include <raptor/argument_parsing/validators.hpp>
#include <raptor/index.hpp>
#include <raptor/update/update.hpp>

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

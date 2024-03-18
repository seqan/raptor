// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::chopper_layout.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <chopper/input_functor.hpp>
#include <chopper/layout/execute.hpp>
#include <chopper/set_up_parser.hpp>
#include <chopper/sketch/check_filenames.hpp>
#include <chopper/sketch/read_data_file.hpp>

#include <raptor/argument_parsing/init_shared_meta.hpp>

namespace raptor
{

void chopper_layout(sharg::parser & parser)
{
    chopper::configuration config;
    set_up_parser(parser, config);
    init_shared_meta(parser);
    parser.info.synopsis.front().insert(0, "raptor layout");
    parser.info.author = "Svenja Mehringer";
    parser.info.email = "svenja.mehringer@fu-berlin.de";

    parser.parse();

    if (!parser.is_option_set("window"))
        config.window_size = config.k;
    else if (config.k > config.window_size)
        throw sharg::parser_error{"The k-mer size cannot be bigger than the window size."};

    config.disable_sketch_output = !parser.is_option_set("output-sketches-to");

    if (std::filesystem::is_empty(config.data_file))
        throw sharg::parser_error{"The input file is empty."};

    std::vector<std::vector<std::string>> filenames{};
    chopper::sketch::read_data_file(config, filenames);

    chopper::sketch::check_filenames(filenames, config);

    config.hibf_config.input_fn = chopper::input_functor{.filenames = filenames,
                                                         .input_are_precomputed_files = config.precomputed_files,
                                                         .kmer_size = config.k,
                                                         .window_size = config.window_size};
    config.hibf_config.number_of_user_bins = filenames.size();

    chopper::layout::execute(config, filenames);
}

} // namespace raptor

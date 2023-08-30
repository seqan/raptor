// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

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
    parser.info.author = "Svenja Mehringer";
    parser.info.email = "svenja.mehringer@fu-berlin.de";

    parser.parse();
    config.disable_sketch_output = !parser.is_option_set("output-sketches-to");

    if (std::filesystem::is_empty(config.data_file))
        throw sharg::parser_error{"The input file is empty."};

    std::vector<std::string> filenames{};
    chopper::sketch::read_data_file(config, filenames);

    chopper::sketch::check_filenames(filenames, config);

    config.hibf_config.input_fn = chopper::input_functor{filenames, config.precomputed_files, config.k};
    config.hibf_config.number_of_user_bins = filenames.size();

    chopper::layout::execute(config, filenames);
}

} // namespace raptor

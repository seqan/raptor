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

#include <chopper/configuration.hpp>
#include <chopper/count/execute.hpp>
#include <chopper/data_store.hpp>
#include <chopper/layout/execute.hpp>
#include <chopper/set_up_parser.hpp>
#include <chopper/sketch/estimate_kmer_counts.hpp>

#include <raptor/argument_parsing/init_shared_meta.hpp>
#include <raptor/layout/raptor_layout.hpp>

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

    // The output streams facilitate writing the layout file in hierarchical structure.
    // chopper::layout::execute currently writes the filled buffers to the output file.
    std::stringstream output_buffer;
    std::stringstream header_buffer;

    chopper::data_store store{.false_positive_rate = config.false_positive_rate,
                              .output_buffer = &output_buffer,
                              .header_buffer = &header_buffer};

    chopper::count::execute(config, store);
    chopper::sketch::estimate_kmer_counts(store);
    chopper::layout::execute(config, store);
}

} // namespace raptor

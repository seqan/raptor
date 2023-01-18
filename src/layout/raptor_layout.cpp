// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <chopper/count/execute.hpp>
#include <chopper/detail_apply_prefix.hpp>
#include <chopper/layout/execute.hpp>
#include <chopper/set_up_parser.hpp>

#include <raptor/layout/raptor_layout.hpp>

namespace raptor
{

void chopper_layout(sharg::parser & parser)
{
    chopper::configuration config;
    set_up_parser(parser, config);

    parser.parse();

    config.input_prefix = config.output_prefix;

    chopper::detail::apply_prefix(config.output_prefix, config.count_filename, config.sketch_directory);

    chopper::count::execute(config);
    chopper::layout::execute(config);
}

} // namespace raptor

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

    try
    {
        parser.parse();
    }
    catch (sharg::parser_error const & ext) // the user did something wrong
    {
        std::cerr << "[CHOPPER ERROR] " << ext.what() << '\n'; // customize your error message
        return;
    }

    config.input_prefix = config.output_prefix;

    chopper::detail::apply_prefix(config.output_prefix, config.count_filename, config.sketch_directory);

    int exit_code{};

    try
    {
        exit_code |= chopper::count::execute(config);
        exit_code |= chopper::layout::execute(config);
    }
    catch (sharg::parser_error const & ext)
    {
        std::cerr << "[CHOPPER ERROR] " << ext.what() << '\n';
    }
}

} // namespace raptor

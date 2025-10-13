// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::chopper_layout.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <string> // for basic_string
#include <vector> // for vector

#include <sharg/auxiliary.hpp> // for parser_meta_data
#include <sharg/parser.hpp>    // for parser

#include <chopper/chopper_layout.hpp> // for chopper_layout
#include <chopper/configuration.hpp>  // for configuration
#include <chopper/set_up_parser.hpp>  // for set_up_parser

namespace raptor
{

void chopper_layout(sharg::parser & parser)
{
    chopper::configuration config;
    set_up_parser(parser, config);
    parser.info.synopsis.front().insert(0, "raptor layout");
    parser.info.author = "Svenja Mehringer";
    parser.info.email = "svenja.mehringer@fu-berlin.de";

    chopper::chopper_layout(config, parser);
}

} // namespace raptor

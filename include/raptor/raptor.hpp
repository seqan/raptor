// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor callable.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <sharg/parser.hpp>

namespace raptor
{

void init_top_level_parser(sharg::parser & parser);
void run_build(sharg::parser & parser);
void run_search(sharg::parser & parser);

} // namespace raptor

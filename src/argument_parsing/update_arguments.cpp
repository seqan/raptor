// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::update_arguments.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <iostream>

#include <raptor/argument_parsing/update_arguments.hpp>

namespace raptor
{

void update_arguments::print_stats() const
{
    std::cerr << "Inserted User Bins: " << user_bins_to_insert.size() << '\n';
    std::cerr << "Full Rebuilds: " << full_rebuilds << '\n';
    std::cerr << "Partial Rebuilds: " << partial_rebuilds << '\n';
}

} // namespace raptor

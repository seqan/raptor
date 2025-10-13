// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::compute_bin_size and raptor::max_bin_count.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <cstddef> // for size_t

#include <raptor/argument_parsing/build_arguments.hpp>   // for build_arguments
#include <raptor/argument_parsing/upgrade_arguments.hpp> // for upgrade_arguments

namespace raptor
{

size_t compute_bin_size(raptor::build_arguments const & arguments);
size_t max_bin_count(upgrade_arguments const & arguments);

} // namespace raptor

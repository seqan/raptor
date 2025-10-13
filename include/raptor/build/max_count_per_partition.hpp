// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::max_count_per_partition.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <cstddef> // for size_t
#include <vector>  // for vector

#include <raptor/argument_parsing/build_arguments.hpp>   // for build_arguments
#include <raptor/argument_parsing/upgrade_arguments.hpp> // for upgrade_arguments
#include <raptor/build/partition_config.hpp>             // for partition_config

namespace raptor
{

std::vector<size_t> max_count_per_partition(partition_config const & cfg, build_arguments const & arguments);
std::vector<size_t> max_count_per_partition(partition_config const & cfg, upgrade_arguments const & arguments);

} // namespace raptor

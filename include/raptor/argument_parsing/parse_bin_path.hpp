// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::parse_bin_path.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <filesystem> // for path
#include <string>     // for string
#include <vector>     // for vector

#include <raptor/argument_parsing/build_arguments.hpp>   // for build_arguments
#include <raptor/argument_parsing/prepare_arguments.hpp> // for prepare_arguments
#include <raptor/argument_parsing/upgrade_arguments.hpp> // for upgrade_arguments

namespace raptor
{

namespace detail
{

void parse_bin_path(std::filesystem::path const & bin_file,
                    std::vector<std::vector<std::string>> & bin_path,
                    bool const is_hibf);

} // namespace detail

void parse_bin_path(build_arguments & arguments);
void parse_bin_path(prepare_arguments & arguments);
void parse_bin_path(upgrade_arguments & arguments);

} // namespace raptor

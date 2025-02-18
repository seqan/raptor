// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::insert_user_bin.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <raptor/argument_parsing/update_arguments.hpp>
#include <raptor/index.hpp>

namespace raptor
{

void insert_user_bin(update_arguments const & arguments, raptor_index<index_structure::hibf> & index);

} // namespace raptor

// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::raptor_update.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <raptor/argument_parsing/update_arguments.hpp>

namespace raptor
{

/*!\brief Updates an existing HIBF index and stores the result.
 * \param[in] arguments The update arguments.
 * \details
 * Loads the index from `arguments.index_file`, applies `arguments.user_bins_to_delete` and
 * `arguments.user_bins_to_insert` in that order, and stores the result at `arguments.out_path`. The input index is
 * left untouched.
 */
void raptor_update(update_arguments const & arguments);

} // namespace raptor

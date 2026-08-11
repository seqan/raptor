// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
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

/*!\brief Inserts the user bins listed in `arguments.user_bins_to_insert` into `index`.
 * \param[in] arguments The update arguments.
 * \param[in,out] index The index to insert into.
 * \details
 * Each user bin is inserted individually. For each one, a suitable IBF and a range of technical bins is determined
 * (see raptor::detail::get_location), the k-mers are inserted there and into all merged bins on the path to the
 * top-level IBF, and the index is rebuilt if that insertion degraded it too much:
 *
 * * If an IBF holds more technical bins than `tmax` allows, it is rebuilt. For the top-level IBF this means a full
 *   rebuild, for any other IBF a partial rebuild of the subtree it roots.
 * * Otherwise, if inserting exceeded the target FPR of a bin, the subtree containing that bin is rebuilt partially,
 *   or fully if the bin lies in the top-level IBF.
 *
 * A partial rebuild recomputes the layout of one subtree, a full rebuild recomputes the whole index from all user
 * bins, including the ones not inserted yet. Consequently, a full rebuild finishes the insertion of every remaining
 * user bin, and this function returns right after it. Setting `arguments.disable_partial_rebuild` replaces every
 * partial rebuild by a full one.
 */
void insert_user_bin(update_arguments const & arguments, raptor_index<index_structure::hibf> & index);

} // namespace raptor

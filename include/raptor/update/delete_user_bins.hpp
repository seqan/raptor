// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::delete_user_bins.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <raptor/argument_parsing/update_arguments.hpp>
#include <raptor/index.hpp>

namespace raptor
{

/*!\brief Deletes the user bins listed in `arguments.user_bins_to_delete` from `index`.
 * \param[in] arguments The update arguments.
 * \param[in,out] index The index to delete from.
 * \details
 * The technical bins holding a deleted user bin are cleared and marked as seqan::hibf::bin_kind::deleted, and their
 * occupancy is reset. The freed technical bins are reused by a subsequent raptor::insert_user_bin.
 *
 * An IBF whose bins all became empty is additionally cleared in its parent's merged bin, and the subtree it roots
 * is retired, see raptor::detail::retire_ibf: their bit vectors are released, but their slots in `ibf_vector` are
 * kept so that the IBF indices stored throughout the index stay valid.
 *
 * User bin IDs are not reassigned, hence the IDs of the remaining user bins stay valid.
 */
void delete_user_bins(update_arguments const & arguments, raptor_index<index_structure::hibf> & index);

} // namespace raptor

// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides helpers to traverse and retire subtrees of an HIBF.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include <raptor/index.hpp>

namespace raptor::detail
{

//!\brief Collects the indices of all IBFs in the subtree rooted at `ibf_idx` into `result`, `ibf_idx` first.
inline void
collect_subtree(std::vector<size_t> & result, raptor_index<index_structure::hibf> const & index, size_t const ibf_idx)
{
    result.push_back(ibf_idx);
    auto const & user_bin_ids = index.ibf().ibf_bin_to_user_bin_id[ibf_idx];
    for (size_t i = 0; i < user_bin_ids.size(); ++i)
    {
        if (user_bin_ids[i] == seqan::hibf::bin_kind::merged)
            collect_subtree(result, index, index.ibf().next_ibf_id[ibf_idx][i]);
    }
}

/*!\brief Retires the IBF at `ibf_idx`: releases its data and marks it as deleted.
 * \details
 * All data of the IBF is released. Its slot in `ibf_vector` is kept, because erasing it would shift, and thereby
 * invalidate, every larger IBF index stored in the index. A retired IBF is recognised by its `prev_ibf_id`, see
 * raptor::detail::is_deleted, and is skipped when looking for a place to insert a user bin.
 *
 * The caller has to unmark the merged bin the IBF used to hang in, and to restore that bin's `next_ibf_id` entry.
 */
inline void retire_ibf(raptor_index<index_structure::hibf> & index, size_t const ibf_idx)
{
    auto & hibf = index.ibf();
    // Assigning an empty container releases the memory; `clear()` would keep the capacity.
    hibf.ibf_vector[ibf_idx] = seqan::hibf::interleaved_bloom_filter{};
    hibf.next_ibf_id[ibf_idx] = std::vector<uint64_t>{};
    hibf.prev_ibf_id[ibf_idx] = {seqan::hibf::bin_kind::deleted, seqan::hibf::bin_kind::deleted};
    hibf.ibf_bin_to_user_bin_id[ibf_idx] = std::vector<uint64_t>{};
}

} // namespace raptor::detail

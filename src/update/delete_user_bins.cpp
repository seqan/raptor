// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::delete_user_bins.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <raptor/update/delete_user_bins.hpp>

#include "subtree.hpp"

namespace raptor
{

void delete_user_bins(update_arguments const & arguments, raptor_index<index_structure::hibf> & index)
{
    std::vector<seqan::hibf::bin_index> technical_bins_to_delete{};
    auto & ibf_bin_to_user_bin_id = index.ibf().ibf_bin_to_user_bin_id;

    for (size_t ibf_index = 0; ibf_index < ibf_bin_to_user_bin_id.size(); ++ibf_index)
    {
        technical_bins_to_delete.clear();

        auto & ibf_index_to_user_bin_id = ibf_bin_to_user_bin_id[ibf_index];
        auto & ibf = index.ibf().ibf_vector[ibf_index];

        for (size_t bin_index = 0; bin_index < ibf_index_to_user_bin_id.size(); ++bin_index)
        {
            uint64_t & user_bin_id = ibf_index_to_user_bin_id[bin_index];
            // This also ensures that invalid user bin IDs are not processed. TODO: Warning/Check?
#ifdef __cpp_lib_ranges_contains
            if (std::ranges::contains(arguments.user_bins_to_delete, user_bin_id))
#else
            if (std::ranges::find(arguments.user_bins_to_delete, user_bin_id)
                != std::ranges::end(arguments.user_bins_to_delete))
#endif
            {
                technical_bins_to_delete.push_back({bin_index});
                user_bin_id = seqan::hibf::bin_kind::deleted;
            }
        }

        if (!technical_bins_to_delete.empty())
        {
            ibf.clear(technical_bins_to_delete);
            for (auto const technical_bin_index : technical_bins_to_delete)
            {
                ibf.occupancy[technical_bin_index.value] = 0u;
            }
            bool const all_zero = std::ranges::all_of(ibf.occupancy,
                                                      [](size_t value)
                                                      {
                                                          return value == 0u;
                                                      });
            // Delete in parent
            if (ibf_index != 0 && all_zero)
            {
                auto const parent = index.ibf().prev_ibf_id[ibf_index];

                auto & parent_ibf = index.ibf().ibf_vector[parent.ibf_idx];
                parent_ibf.clear(seqan::hibf::bin_index{parent.bin_idx});

                parent_ibf.occupancy[parent.bin_idx] = 0u;
                ibf_bin_to_user_bin_id[parent.ibf_idx][parent.bin_idx] = seqan::hibf::bin_kind::deleted;
                // The bin is no longer a merged bin. Every reader checks `ibf_bin_to_user_bin_id` before following
                // `next_ibf_id`, but seqan::hibf documents `next_ibf_id[i][b] == i` for a bin that is not merged.
                index.ibf().next_ibf_id[parent.ibf_idx][parent.bin_idx] = parent.ibf_idx;

                // The subtree is no longer reachable: its merged bin in the parent is gone. Retire it, otherwise
                // raptor::detail::max_ibf_sizes keeps offering these IBFs as a place to insert into, and a rebuild
                // triggered by such an insertion would look for a merged bin that no longer exists.
                //
                // Alternative: keep the bin merged instead of deleting it, i.e. drop the assignment above. The
                // subtree then stays reachable but empty; a query never descends into it because the merged bin was
                // cleared, and inserting into it repopulates the merged bin. That keeps the empty IBFs around, but
                // makes them reusable as free space.
                //
                // TODO: Retiring a subtree can leave its parent empty in turn. This loop only moves forward, so a
                //       parent with a smaller index is not revisited. Such an IBF stays alive but empty, which
                //       costs memory but is otherwise harmless: it is a valid place to insert into.
                std::vector<size_t> subtree{};
                detail::collect_subtree(subtree, index, ibf_index);
                for (size_t const subtree_ibf_index : subtree)
                    detail::retire_ibf(index, subtree_ibf_index);
            }
        }
    }
}

} // namespace raptor

// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::delete_user_bins.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <algorithm> // for __all_of, __contains, all_of, contains
#include <cstddef>   // for size_t
#include <cstdint>   // for uint64_t
#include <iostream>  // for basic_ostream, operator<<, cerr, basic_ios
#include <vector>    // for vector
#include <version>   // for __cpp_lib_ranges_contains

#include <hibf/hierarchical_interleaved_bloom_filter.hpp> // for hierarchical_interleaved_bloom_filter, deleted
#include <hibf/interleaved_bloom_filter.hpp>              // for bin_index, interleaved_bloom_filter

#include <raptor/argument_parsing/update_arguments.hpp> // for update_arguments
#include <raptor/index.hpp>                             // for raptor_index, hibf
#include <raptor/update/delete_user_bins.hpp>           // for delete_user_bins

namespace raptor
{

void delete_user_bins(update_arguments const & arguments, raptor_index<index_structure::hibf> & index)
{
    std::cerr << "\nDeleting user bins: ";
    for (auto const bin : arguments.user_bins_to_delete)
        std::cerr << bin << ' ';
    std::cerr << '\n';

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
            }
        }
    }
}

} // namespace raptor

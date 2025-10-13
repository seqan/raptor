// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::insert_user_bin.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <cstddef> // for size_t
#include <cstdint> // for uint64_t
#include <limits>  // for numeric_limits
#include <vector>  // for vector

#include <hibf/build/insert_into_ibf.hpp>                 // for insert_into_ibf
#include <hibf/contrib/robin_hood.hpp>                    // for unordered_flat_set
#include <hibf/hierarchical_interleaved_bloom_filter.hpp> // for hierarchical_interleaved_bloom_filter
#include <hibf/interleaved_bloom_filter.hpp>              // for interleaved_bloom_filter

#include <raptor/index.hpp> // for raptor_index, hibf

#include "is_fpr_exceeded.hpp" // for is_fpr_exceeded
#include "strong_types.hpp"    // for insert_location, rebuild_location

namespace raptor::detail
{

void insert_into_ibf(robin_hood::unordered_flat_set<uint64_t> const & kmers,
                     insert_location const & insert_location,
                     raptor_index<index_structure::hibf> & index,
                     rebuild_location & rebuild_index_tuple)
{
    auto & ibf = index.ibf().ibf_vector[insert_location.ibf_idx];

    // std::cout << "\nIBF " << insert_location.ibf_idx << '\n';
    // std::cout << "Bin " << insert_location.bin_idx << '\n';
    // std::cout << "Insert " << insert_location.number_of_bins << " bins\n";

    // std::cout << "Kmer content " << kmers.size() << '\n';
    // std::cout << "FPR before " << compute_fpr(ibf, insert_location.bin_idx) << '\n';
    // std::cout << "Occupancy before " << ibf.occupancy[insert_location.bin_idx] << '\n';

    seqan::hibf::build::insert_into_ibf(kmers,
                                        insert_location.number_of_bins,
                                        insert_location.bin_idx,
                                        ibf,
                                        index.ibf().fill_ibf_timer);

    // std::cout << "Occupancy after " << ibf.occupancy[insert_location.bin_idx] << '\n';

    // TODO: Won't the kmers be evenly split? In this case, one computation is enough.
    // for (size_t i = insert_location.bin_idx; i < insert_location.bin_idx + insert_location.number_of_bins; ++i)
    // {
    if (is_fpr_exceeded(index, insert_location))
    {
        rebuild_index_tuple.ibf_idx = insert_location.ibf_idx;
        rebuild_index_tuple.bin_idx = insert_location.bin_idx;
    }
    // }

    // TODO
    // // to improve the implementation, Perhaps do the FPR calculations for all bins to which kmers will be inserted before actually inserting.
    // index.ibf().update_occupancy_table(kmers.size()-union_count, ibf_idx, start_bin_idx, number_of_bins);
    // auto fpr = index.ibf().update_fpr(ibf_idx, start_bin_idx, number_of_bins); // this should be done after updating the occupancy table.

    // // Keep: If any fpr is too high, the index needs to be rebuild.
    // if (fpr > index.ibf().fpr_max){
    //     //assert(index.ibf().next_ibf_id[ibf_idx][start_bin_idx] != ibf_idx); //assert that it is not a leaf bin. fpr should not reach fpr max for the leaf bin ibf, because we search a location such that it is 'feasible', i.e. binsize should be sufficient to accomodate new UB. However, this is not the case if one also does sequence insertions.
    //     rebuild_index_tuple = std::make_tuple(ibf_idx, start_bin_idx);
    // }
}

rebuild_location insert_tb_and_parents(robin_hood::unordered_flat_set<uint64_t> const & kmers,
                                       insert_location insert_location,
                                       raptor_index<index_structure::hibf> & index)
{
    rebuild_location rebuild_location{.ibf_idx = std::numeric_limits<size_t>::max(),
                                      .bin_idx = std::numeric_limits<size_t>::max()};
    while (true)
    {
        insert_into_ibf(kmers, insert_location, index, rebuild_location);
        if (insert_location.ibf_idx == 0u)
            break;
        auto const parent = index.ibf().prev_ibf_id[insert_location.ibf_idx];
        insert_location.ibf_idx = parent.ibf_idx;
        insert_location.bin_idx = parent.bin_idx;
        insert_location.number_of_bins = 1u; // merged bin
    }

    return rebuild_location;
}

} // namespace raptor::detail

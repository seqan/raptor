// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::detail::insert_tb_and_parents.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <optional>

#include <hibf/build/insert_into_ibf.hpp>

#include <raptor/index.hpp>

#include "is_fpr_exceeded.hpp"
#include "strong_types.hpp"

namespace raptor::detail
{

/*!\brief Inserts `kmers` into the technical bins given by `insert_location`.
 * \param[in] kmers The k-mers to insert.
 * \param[in] insert_location Where to insert.
 * \param[in,out] index The index to insert into.
 * \param[in,out] rebuild_index_tuple Set to `insert_location` if the insertion pushed the bin past its target FPR.
 * \details The k-mers are distributed evenly if `insert_location` spans more than one technical bin.
 */
void insert_into_ibf(robin_hood::unordered_flat_set<uint64_t> const & kmers,
                     insert_location const & insert_location,
                     raptor_index<index_structure::hibf> & index,
                     std::optional<rebuild_location> & rebuild_index_tuple)
{
    auto & ibf = index.ibf().ibf_vector[insert_location.ibf_idx];

    seqan::hibf::build::insert_into_ibf(kmers,
                                        insert_location.number_of_bins,
                                        insert_location.bin_idx,
                                        ibf,
                                        index.ibf().fill_ibf_timer);

    if (is_fpr_exceeded(index, insert_location))
        rebuild_index_tuple = rebuild_location{.ibf_idx = insert_location.ibf_idx, .bin_idx = insert_location.bin_idx};
}

/*!\brief Inserts `kmers` at `insert_location` and into every merged bin on the path to the top-level IBF.
 * \param[in] kmers The k-mers to insert.
 * \param[in] insert_location Where the user bin is stored.
 * \param[in,out] index The index to insert into.
 * \returns The bin closest to the top-level IBF whose target FPR is exceeded, or `std::nullopt` if none is.
 * \details
 * A query only reaches a lower-level IBF through the merged bins above it, hence the k-mers have to be present in
 * all of them. Walking upwards means a returned location is the highest one that exceeds its FPR, i.e. rebuilding
 * there also covers every exceeded bin below it.
 */
std::optional<rebuild_location> insert_tb_and_parents(robin_hood::unordered_flat_set<uint64_t> const & kmers,
                                                      insert_location insert_location,
                                                      raptor_index<index_structure::hibf> & index)
{
    std::optional<rebuild_location> rebuild_location{};
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

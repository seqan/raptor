// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::detail::is_fpr_exceeded.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>

#include <raptor/index.hpp>

#include "strong_types.hpp"

namespace raptor::detail
{

//!\brief The level of the IBF at `ibf_idx` within the HIBF. The top-level IBF has level 0.
inline size_t level_of(raptor_index<index_structure::hibf> const & index, size_t ibf_idx)
{
    [[maybe_unused]] size_t const number_of_ibfs = index.ibf().ibf_vector.size();
    size_t level{};

    while (ibf_idx != 0u)
    {
        assert(level < number_of_ibfs && "Cycle in the chain of parent IBFs.");
        ++level;
        ibf_idx = index.ibf().prev_ibf_id[ibf_idx].ibf_idx;
    }

    return level;
}

/*!\brief Whether the technical bin at `bin_idx` in the IBF at `ibf_idx` exceeds its target FPR.
 * \param[in] index The index.
 * \param[in] ibf_idx The IBF containing the bin.
 * \param[in] bin_idx The bin to check.
 * \param[in] level The level of `ibf_idx`; only relevant for merged bins.
 */
inline bool is_fpr_exceeded_impl(raptor_index<index_structure::hibf> const & index,
                                 size_t const ibf_idx,
                                 size_t const bin_idx,
                                 size_t const level)
{
    auto const & hibf = index.ibf();
    auto const & ibf = hibf.ibf_vector[ibf_idx];

    double const new_fpr = [&ibf, bin_idx]()
    {
        double const exp_arg =
            (ibf.hash_function_count() * ibf.occupancy[bin_idx]) / static_cast<double>(ibf.bin_size());
        double const log_arg = 1.0 - std::exp(-exp_arg);
        return std::exp(ibf.hash_function_count() * std::log(log_arg));
    }();

    bool const is_bin_merged = hibf.ibf_bin_to_user_bin_id[ibf_idx][bin_idx] == seqan::hibf::bin_kind::merged;
    double const target_fpr = [&]() -> double
    {
        // A bin holding a user bin must respect the maximum FPR.
        if (!is_bin_merged)
            return index.fpr();

        // A merged bin only needs to respect the relaxed FPR, because a query passing it is checked again in the
        // IBF below. Rebuilding is more expensive the closer an IBF is to the top, so merged bins are granted
        // additional slack that decays with the level: 2x the relaxed FPR at the top level, 1.25x one level below,
        // and so on. The target FPR never exceeds 0.95, or the relaxed FPR itself if that is already higher.
        double const relaxed_fpr = index.config().relaxed_fpr;
        double const factor = 1.0 + 1.0 / ((level + 1u) * (level + 1u));
        return std::min(relaxed_fpr * factor, std::max(relaxed_fpr, 0.95));
    }();

    return new_fpr > target_fpr;
}

//!\brief Whether the first bin of `insert_location` exceeds its target FPR.
inline bool is_fpr_exceeded(raptor_index<index_structure::hibf> const & index, insert_location const & insert_location)
{
    return is_fpr_exceeded_impl(index,
                                insert_location.ibf_idx,
                                insert_location.bin_idx,
                                level_of(index, insert_location.ibf_idx));
}

//!\brief Whether the bin of `rebuild_location` exceeds its target FPR.
inline bool is_fpr_exceeded(raptor_index<index_structure::hibf> const & index,
                            rebuild_location const & rebuild_location)
{
    auto const level = level_of(index, rebuild_location.ibf_idx);
    assert(level == 0u); // This overload should only be called with the top-level IBF as rebuild location.
    return is_fpr_exceeded_impl(index, rebuild_location.ibf_idx, rebuild_location.bin_idx, level);
}

} // namespace raptor::detail

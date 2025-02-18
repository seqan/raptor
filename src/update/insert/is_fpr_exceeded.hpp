// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::insert_user_bin.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <cmath>
#include <cstddef>

#include <raptor/index.hpp>

#include "strong_types.hpp"

namespace raptor::detail
{

inline bool is_fpr_exceeded_impl(raptor_index<index_structure::hibf> const & index,
                                 size_t const ibf_idx,
                                 size_t const bin_idx,
                                 bool const is_toplevel)
{
    auto & hibf = index.ibf();
    auto & ibf = hibf.ibf_vector[ibf_idx];

    double const new_fpr = [&ibf, &bin_idx]()
    {
        double const exp_arg =
            (ibf.hash_function_count() * ibf.occupancy[bin_idx]) / static_cast<double>(ibf.bin_size());
        double const log_arg = 1.0 - std::exp(-exp_arg);
        return std::exp(ibf.hash_function_count() * std::log(log_arg));
    }();

    bool const is_bin_merged = hibf.ibf_bin_to_user_bin_id[ibf_idx][bin_idx] == seqan::hibf::bin_kind::merged;
    double const target_fpr = [&]() -> double
    {
        if (!is_bin_merged)
            return index.fpr();

        double const relaxed_fpr = index.config().relaxed_fpr;
        return relaxed_fpr * (is_toplevel ? std::min(relaxed_fpr * 1.25, std::max(relaxed_fpr, 0.95)) : 1.0);
    }();

    return new_fpr > target_fpr;
}

inline bool is_fpr_exceeded(raptor_index<index_structure::hibf> const & index, insert_location const & insert_location)
{
    return is_fpr_exceeded_impl(index, insert_location.ibf_idx, insert_location.bin_idx, false);
}

inline bool is_fpr_exceeded(raptor_index<index_structure::hibf> const & index,
                            rebuild_location const & rebuild_location)
{
    return is_fpr_exceeded_impl(index, rebuild_location.ibf_idx, rebuild_location.bin_idx, true);
}

} // namespace raptor::detail

// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::insert_user_bin.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <cstddef>

namespace raptor::detail
{

struct bookkeeping_arguments
{
    size_t ibf_idx;
    size_t old_number_of_bins;
    size_t number_of_new_bins;
};

struct insert_location
{
    size_t ibf_idx;
    size_t bin_idx;
    size_t number_of_bins;
};

struct rebuild_location
{
    size_t ibf_idx;
    size_t bin_idx;
};

struct max_elements_parameters
{
    double fpr{};
    size_t hash_count{};
    size_t bin_size{};
};

struct ibf_max
{
    size_t max_elements;
    size_t ibf_idx;

    constexpr auto operator<=>(ibf_max const & other) const = default;
};

struct required_technical_bins_parameters
{
    size_t bin_size{};
    size_t elements{};
    double fpr{};
    size_t hash_count{};
    size_t max_elements{};
};

struct ibf_location
{
    size_t ibf_idx;
    size_t max_elements;
};

} // namespace raptor::detail

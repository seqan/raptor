// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides the helper types shared by the raptor::insert_user_bin implementation.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <cstddef>

namespace raptor::detail
{

/*!\brief The number of technical bins a lower-level IBF may exceed `tmax` by before it is rebuilt.
 * \details Without this slack, an IBF that reached `tmax` would trigger a rebuild on every subsequent insertion.
 */
inline constexpr size_t tmax_slack{64u};

//!\brief Arguments for raptor::detail::update_bookkeeping.
struct bookkeeping_arguments
{
    //!\brief The IBF the new user bin was placed in.
    size_t ibf_idx;
    //!\brief The first technical bin that the new user bin occupies.
    size_t first_bin_idx;
    //!\brief The number of technical bins the new user bin occupies.
    size_t number_of_new_bins;
};

//!\brief A range of consecutive technical bins within one IBF, i.e. where a user bin is stored.
struct insert_location
{
    size_t ibf_idx;
    //!\brief The first technical bin of the range.
    size_t bin_idx;
    //!\brief The number of technical bins in the range. Greater than 1 if the user bin is split.
    size_t number_of_bins;
};

//!\brief A merged bin whose subtree needs to be rebuilt.
struct rebuild_location
{
    //!\brief The IBF containing the merged bin, i.e. the parent of the subtree to rebuild.
    size_t ibf_idx;
    //!\brief The merged bin within that IBF.
    size_t bin_idx;
};

//!\brief Arguments for raptor::detail::max_elements.
struct max_elements_parameters
{
    //!\brief The false positive rate that must not be exceeded.
    double fpr{};
    //!\brief The number of hash functions of the IBF.
    size_t hash_count{};
    //!\brief The size of a single technical bin in bits.
    size_t bin_size{};
};

/*!\brief The maximum number of elements a technical bin of an IBF can hold.
 * \details Sorts by `max_elements` first, hence a sorted range of these can be binary searched by k-mer count.
 */
struct ibf_max
{
    //!\brief The number of elements a technical bin of this IBF holds at exactly the maximum FPR.
    size_t max_elements;
    //!\brief The IBF this refers to.
    size_t ibf_idx;

    constexpr auto operator<=>(ibf_max const & other) const = default;
};

//!\brief Arguments for raptor::detail::required_technical_bins.
struct required_technical_bins_parameters
{
    //!\brief The size of a single technical bin in bits.
    size_t bin_size{};
    //!\brief The number of elements to store.
    size_t elements{};
    //!\brief The false positive rate that must not be exceeded.
    double fpr{};
    //!\brief The number of hash functions of the IBF.
    size_t hash_count{};
    //!\brief The number of elements a single technical bin of this IBF can hold.
    size_t max_elements{};
};

} // namespace raptor::detail

// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::detail::get_location.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <optional>
#include <vector>

#include <hibf/misc/divide_and_ceil.hpp>

#include <raptor/index.hpp>

#include "strong_types.hpp"

namespace raptor::detail
{

/*!\brief The number of technical bins needed to store `params.elements` elements without exceeding `params.fpr`.
 * \details
 * Splitting a user bin over several technical bins lowers the load of each one, but a query has to pass all of
 * them, so the combined FPR shrinks more slowly than the individual one. Hence the number of bins is increased
 * until the combined FPR is within budget.
 *
 * Mirrors `seqan::hibf::hierarchical_interleaved_bloom_filter::number_of_bins`.
 */
size_t required_technical_bins(required_technical_bins_parameters const & params)
{
    if (params.elements <= params.max_elements)
        return 1uz;

    auto compute_fpr = [&](size_t const elements)
    {
        double const exp_arg = (params.hash_count * elements) / static_cast<double>(params.bin_size);
        double const log_arg = 1.0 - std::exp(-exp_arg);
        return std::exp(params.hash_count * std::log(log_arg));
    };

    auto compute_split_fpr = [&](size_t const split)
    {
        double const fpr_tb = compute_fpr(seqan::hibf::divide_and_ceil(params.elements, split));
        return 1.0 - std::exp(std::log1p(-fpr_tb) * split);
    };

    size_t number_of_bins = seqan::hibf::divide_and_ceil(params.elements, params.max_elements);

    while (compute_split_fpr(number_of_bins) > params.fpr)
        ++number_of_bins;

    return number_of_bins;
}

//!\brief How far an IBF may be grown while looking for space for a new user bin.
enum class growth_policy : uint8_t
{
    //!\brief Only use space that is already allocated.
    none,
    //!\brief Grow, but at most to `tmax + tmax_slack` technical bins.
    limited,
    //!\brief Grow without a limit. The IBF is marked as resized, hence it is grown at most once.
    unlimited
};

/*!\brief Claims `number_of_bins` consecutive technical bins in the IBF at `ibf_idx`.
 * \returns The index of the first claimed bin, or `std::nullopt` if `policy` does not permit claiming any.
 * \details
 * Growing the IBF is a side effect of claiming bins; `policy` decides whether that is permitted. The top-level IBF
 * is never grown here, because doing so requires a full rebuild; raptor::detail::find_insert_location handles it.
 */
std::optional<size_t> try_claim_bins(raptor_index<index_structure::hibf> & index,
                                     size_t const ibf_idx,
                                     size_t const number_of_bins,
                                     growth_policy const policy)
{
    auto & ibf = index.ibf().ibf_vector[ibf_idx];

    // Without occupancy tracking, `occupancy` is all-zero and an empty run says nothing about the bins actually
    // being empty. raptor::update_parsing rejects such indices up front.
    assert(ibf.track_occupancy);

    // Prefer a run of empty bins, e.g. bins that became empty through deletions. If there is none, the new bins are
    // appended at the end. `occupancy` covers all technical bins, so a run may also lie in the yet unused ones.
    size_t const bin_idx = [&]() -> size_t
    {
        auto const empty_run = std::ranges::search_n(ibf.occupancy, number_of_bins, 0u);
        if (!empty_run.empty())
            return static_cast<size_t>(std::ranges::distance(ibf.occupancy.begin(), empty_run.begin()));
        return ibf.bin_count();
    }();

    size_t const new_bin_count{bin_idx + number_of_bins};

    // The run lies within the bins that are already in use. `try_increase_bin_number_to` would reject this request
    // because it would decrease the bin count.
    if (new_bin_count <= ibf.bin_count())
        return bin_idx;

    // The run extends past the bins in use, but still fits into the allocated technical bins.
    if (ibf.try_increase_bin_number_to(seqan::hibf::bin_count{new_bin_count}))
        return bin_idx;

    // Everything below grows the underlying bit vector.
    if (policy == growth_policy::none || ibf_idx == 0u || index.was_resized(ibf_idx))
        return std::nullopt;

    if (policy == growth_policy::limited && new_bin_count > index.config().tmax + tmax_slack)
        return std::nullopt;

    if (policy == growth_policy::unlimited)
        index.mark_resized(ibf_idx);

    ibf.increase_bin_number_to(seqan::hibf::bin_count{new_bin_count});
    return bin_idx;
}

/*!\brief Determines where a user bin with `kmer_count` k-mers should be inserted.
 * \details
 * Candidates are visited in order of increasing difference between their maximum number of elements and
 * `kmer_count`, first without growing any IBF, then allowing growth up to `tmax + tmax_slack`, and finally allowing
 * unrestricted growth. If no IBF can take the user bin, the top-level IBF is grown.
 * \returns Where to insert the user bin, or `std::nullopt` if not even the top-level IBF can be grown to take it
 *          without exceeding `tmax`, i.e. if a full rebuild is unavoidable.
 */
std::optional<insert_location> find_insert_location(std::vector<ibf_max> const & max_ibf_sizes,
                                                    size_t const kmer_count,
                                                    raptor_index<index_structure::hibf> & index)
{
    assert(!max_ibf_sizes.empty());
    assert(std::ranges::is_sorted(max_ibf_sizes));

    auto try_candidate = [&](ibf_max const & candidate, growth_policy const policy) -> std::optional<insert_location>
    {
        auto const & ibf = index.ibf().ibf_vector[candidate.ibf_idx];
        size_t const number_of_bins = required_technical_bins({.bin_size = ibf.bin_size(),
                                                               .elements = kmer_count,
                                                               .fpr = index.fpr(),
                                                               .hash_count = ibf.hash_function_count(),
                                                               .max_elements = candidate.max_elements});

        if (auto const bin_idx = try_claim_bins(index, candidate.ibf_idx, number_of_bins, policy); bin_idx.has_value())
            return insert_location{.ibf_idx = candidate.ibf_idx,
                                   .bin_idx = bin_idx.value(),
                                   .number_of_bins = number_of_bins};

        return std::nullopt;
    };

    // `max_ibf_sizes` is sorted by `max_elements`, hence visiting candidates in order of increasing
    // `|max_elements - kmer_count|` is an outward expansion around the first IBF that can hold `kmer_count`.
    size_t const split = static_cast<size_t>(
        std::ranges::distance(max_ibf_sizes.begin(),
                              std::ranges::lower_bound(max_ibf_sizes, kmer_count, {}, &ibf_max::max_elements)));

    auto search = [&](growth_policy const policy) -> std::optional<insert_location>
    {
        size_t too_small{split};    // Candidates in [0, too_small) have not been visited yet.
        size_t large_enough{split}; // Candidates in [large_enough, size) have not been visited yet.

        while (too_small != 0u || large_enough != max_ibf_sizes.size())
        {
            bool take_larger{};
            if (too_small == 0u)
                take_larger = true;
            else if (large_enough == max_ibf_sizes.size())
                take_larger = false;
            else
                take_larger = (max_ibf_sizes[large_enough].max_elements - kmer_count)
                            < (kmer_count - max_ibf_sizes[too_small - 1u].max_elements);

            ibf_max const & candidate = take_larger ? max_ibf_sizes[large_enough++] : max_ibf_sizes[--too_small];

            if (auto const result = try_candidate(candidate, policy); result.has_value())
                return result;
        }

        return std::nullopt;
    };

    for (growth_policy const policy : {growth_policy::none, growth_policy::limited, growth_policy::unlimited})
        if (auto const result = search(policy); result.has_value())
            return result;

    // No IBF has space and none of them may be grown any further: grow the top-level IBF.
    auto const top_level = std::ranges::find(max_ibf_sizes, 0uz, &ibf_max::ibf_idx);
    // `max_ibf_sizes` omits deleted IBFs, but the top-level IBF is never deleted: it is not reachable from any
    // merged bin, and hibf initialises its `prev_ibf_id` to `{.ibf_idx = 0, .bin_idx = 0}`.
    assert(top_level != max_ibf_sizes.end());

    auto & ibf = index.ibf().ibf_vector[0];
    size_t const number_of_bins = required_technical_bins({.bin_size = ibf.bin_size(),
                                                           .elements = kmer_count,
                                                           .fpr = index.fpr(),
                                                           .hash_count = ibf.hash_function_count(),
                                                           .max_elements = top_level->max_elements});
    size_t const bin_idx = ibf.bin_count();

    // Growing the top-level IBF past tmax means the whole index has to be rebuilt. Report that instead of growing,
    // so that the caller can rebuild right away rather than inserting the user bin into a discarded index first.
    if (bin_idx + number_of_bins > index.config().tmax)
        return std::nullopt;

    index.mark_resized(0u);
    ibf.increase_bin_number_to(seqan::hibf::bin_count{bin_idx + number_of_bins});

    return insert_location{.ibf_idx = 0u, .bin_idx = bin_idx, .number_of_bins = number_of_bins};
}

/*!\brief Registers a newly placed user bin in the bookkeeping of its IBF.
 * \details
 * Marks the claimed technical bins as occupied, points them at the new user bin, and grows the per-bin vectors if
 * the user bin was appended at the end. The claimed bins may also be reclaimed ones from the middle of the IBF, in
 * which case the vectors already cover them and only their contents change.
 */
void update_bookkeeping(bookkeeping_arguments const & args, raptor_index<index_structure::hibf> & index)
{
    auto & hibf = index.ibf();
    auto & ibf = hibf.ibf_vector[args.ibf_idx];
    auto & next_ibf_id = hibf.next_ibf_id[args.ibf_idx];
    auto & ibf_bin_to_user_bin_id = hibf.ibf_bin_to_user_bin_id[args.ibf_idx];
    assert(next_ibf_id.size() == ibf_bin_to_user_bin_id.size());

    size_t const last_bin_idx = args.first_bin_idx + args.number_of_new_bins;

    // The claimed bins may lie within the IBF, in which case the bookkeeping already covers them.
    if (last_bin_idx > next_ibf_id.size())
    {
        next_ibf_id.resize(last_bin_idx, args.ibf_idx);
        ibf_bin_to_user_bin_id.resize(last_bin_idx, hibf.number_of_user_bins);
    }

    for (size_t i = args.first_bin_idx; i < last_bin_idx; ++i)
    {
        ibf.occupancy[i] = 1u;
        next_ibf_id[i] = args.ibf_idx; // Not a merged bin.
        ibf_bin_to_user_bin_id[i] = hibf.number_of_user_bins;
    }

    ++hibf.number_of_user_bins;
}

/*!\brief Determines where to store a user bin with `kmer_count` k-mers, and reserves the technical bins for it.
 * \param[in] max_ibf_sizes The capacity of each IBF, see raptor::detail::max_ibf_sizes.
 * \param[in] kmer_count The number of k-mers of the user bin.
 * \param[in,out] index The index to insert into.
 * \returns Where to store the user bin, or `std::nullopt` if a full rebuild is unavoidable, in which case `index`
 *          is left unchanged.
 */
std::optional<insert_location> get_location(std::vector<ibf_max> const & max_ibf_sizes,
                                            size_t const kmer_count,
                                            raptor_index<index_structure::hibf> & index)
{
    std::optional<insert_location> const location = find_insert_location(max_ibf_sizes, kmer_count, index);

    if (location.has_value())
        update_bookkeeping({.ibf_idx = location->ibf_idx,
                            .first_bin_idx = location->bin_idx,
                            .number_of_new_bins = location->number_of_bins},
                           index);

    return location;
}

} // namespace raptor::detail

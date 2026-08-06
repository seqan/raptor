// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::insert_user_bin.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <hibf/misc/divide_and_ceil.hpp>
#include <hibf/misc/iota_vector.hpp>

#include <raptor/index.hpp>

#include "strong_types.hpp"

namespace raptor::detail
{

// raptor::hierarchical_interleaved_bloom_filter::number_of_bins(size_t ibf_idx, size_t kmer_count)
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

enum class extension_mode : uint8_t
{
    none,
    once,
    full
};

std::optional<size_t> find_empty_bin_idx(raptor_index<index_structure::hibf> & index,
                                         size_t const ibf_idx,
                                         size_t const number_of_bins,
                                         extension_mode const mode)
{
    auto & ibf = index.ibf().ibf_vector[ibf_idx];
    size_t const ibf_bin_count = [&]() -> size_t
    {
        auto search_result = std::ranges::search_n(ibf.occupancy, number_of_bins, 0u);
        if (!search_result.empty())
            return std::ranges::distance(ibf.occupancy.begin(), search_result.begin());
        return ibf.bin_count();
    }();

    // If nothing has been returned, no appropriate empty bin has been found and the bin idx will be the size of the IBF,
    // BUG: Deleted bins!
    size_t const new_bin_count{ibf_bin_count + number_of_bins};

    // If we can increase the number of bins without resizing the underlying bitvector
    if (ibf.try_increase_bin_number_to(seqan::hibf::bin_count{new_bin_count}))
        return ibf_bin_count;

    if (mode == extension_mode::none || ibf_idx == 0uz || index.is_resized[ibf_idx])
        return std::nullopt;

    if (mode == extension_mode::once && new_bin_count <= index.config().tmax + 64uz)
    {
        ibf.increase_bin_number_to(seqan::hibf::bin_count{new_bin_count});
        return ibf_bin_count;
    }

    if (mode == extension_mode::full)
    {
        index.is_resized[ibf_idx] = true;
        ibf.increase_bin_number_to(seqan::hibf::bin_count{new_bin_count});
        return ibf_bin_count;
    }

    return std::nullopt;
}

ibf_location find_ibf_size_splitting(std::vector<ibf_max> const & max_ibf_sizes,
                                     size_t const kmer_count,
                                     raptor_index<index_structure::hibf> & index)
{
    // Create indices and sort by absolute difference from kmer_count
    std::vector<size_t> projection = seqan::hibf::iota_vector(max_ibf_sizes.size());

    std::ranges::sort(projection,
                      [&](size_t a, size_t b)
                      {
                          auto diff_a =
                              std::abs(static_cast<std::ptrdiff_t>(max_ibf_sizes[a].max_elements - kmer_count));
                          auto diff_b =
                              std::abs(static_cast<std::ptrdiff_t>(max_ibf_sizes[b].max_elements - kmer_count));
                          return diff_a < diff_b;
                      });

    auto kernel = [&](extension_mode const extend_tmax) -> std::optional<ibf_location>
    {
        for (size_t const projection_idx : projection)
        {
            ibf_max const candidate = max_ibf_sizes[projection_idx];
            auto & ibf = index.ibf().ibf_vector[candidate.ibf_idx];

            size_t const number_of_bins = required_technical_bins({.bin_size = ibf.bin_size(),
                                                                   .elements = kmer_count,
                                                                   .fpr = index.fpr(),
                                                                   .hash_count = ibf.hash_function_count(),
                                                                   .max_elements = candidate.max_elements});

            if (auto const bin_idx = find_empty_bin_idx(index, candidate.ibf_idx, number_of_bins, extend_tmax);
                bin_idx.has_value())
            {
                return ibf_location{.ibf_idx = candidate.ibf_idx,
                                    .bin_idx = bin_idx.value(),
                                    .max_elements = candidate.max_elements};
            }
        }

        return std::nullopt;
    };

    if (auto const result = kernel(extension_mode::none); result.has_value())
        return result.value();

    if (auto const result = kernel(extension_mode::once); result.has_value())
        return result.value();

    if (auto const result = kernel(extension_mode::full); result.has_value())
        return result.value();

    // If there are no empty bins available whatsoever, use the top level IBF.
    auto & top_level_ibf = index.ibf().ibf_vector[0];

    auto top_level_index = [&]()
    {
        auto max_sizes_it = std::ranges::find_if(max_ibf_sizes,
                                                 [](ibf_max const & m)
                                                 {
                                                     return m.ibf_idx == 0;
                                                 });
        assert(max_sizes_it != max_ibf_sizes.end());
        return std::ranges::distance(max_ibf_sizes.begin(), max_sizes_it);
    }();

    size_t const number_of_bins =
        required_technical_bins({.bin_size = top_level_ibf.bin_size(),
                                 .elements = kmer_count,
                                 .fpr = index.fpr(),
                                 .hash_count = top_level_ibf.hash_function_count(),
                                 .max_elements = max_ibf_sizes[top_level_index].max_elements});

    auto bin_idx = find_empty_bin_idx(index, 0, number_of_bins, extension_mode::none);

    return {.ibf_idx = 0,
            .bin_idx = bin_idx.value_or(std::numeric_limits<size_t>::max()),
            .max_elements = max_ibf_sizes[top_level_index].max_elements};
}

void update_bookkeeping(bookkeeping_arguments const & args, raptor_index<index_structure::hibf> & index)
{
    auto & ibf = index.ibf().ibf_vector[args.ibf_idx];
    size_t const new_number_of_bins = args.old_number_of_bins + args.number_of_new_bins;
    for (size_t i = args.old_number_of_bins; i < new_number_of_bins; ++i)
    {
        ibf.occupancy[i] = 1u;
    }
    index.ibf().next_ibf_id[args.ibf_idx].resize(new_number_of_bins, args.ibf_idx);
    index.ibf().ibf_bin_to_user_bin_id[args.ibf_idx].resize(new_number_of_bins, index.ibf().number_of_user_bins);
    index.ibf().number_of_user_bins += 1;
}

insert_location get_location(std::vector<ibf_max> const & max_ibf_sizes,
                             size_t const kmer_count,
                             raptor_index<index_structure::hibf> & index)
{
    auto const [ibf_idx, bidx, max_elements] = find_ibf_size_splitting(max_ibf_sizes, kmer_count, index);
    auto bin_idx = bidx;

    auto & ibf = index.ibf().ibf_vector[ibf_idx];

    size_t const number_of_bins = required_technical_bins({.bin_size = ibf.bin_size(),
                                                           .elements = kmer_count,
                                                           .fpr = index.fpr(),
                                                           .hash_count = ibf.hash_function_count(),
                                                           .max_elements = max_elements});

    // uint64_t bin_idx = find_empty_bin_idx(index, ibf_idx, number_of_bins);

    // TODO: empty bins percentage
    // The current solution resizes the IBF here, but it would be more efficient to check the tmax function after calculating the new bin count and perhaps trigger a partial rebuild dirctly.
    if (bin_idx == std::numeric_limits<size_t>::max())
    {
        index.is_resized[ibf_idx] = true;
        bin_idx = ibf.bin_count();
        // std::cerr << "[DEBUG] Forced increase[" << ibf_idx << "]: " << bin_idx << " to " << bin_idx + number_of_bins << '\n';
        ibf.increase_bin_number_to(seqan::hibf::bin_count{bin_idx + number_of_bins});
    }

    update_bookkeeping({.ibf_idx = ibf_idx, .old_number_of_bins = bin_idx, .number_of_new_bins = number_of_bins},
                       index);

    return insert_location{.ibf_idx = ibf_idx, .bin_idx = bin_idx, .number_of_bins = number_of_bins};
}

} // namespace raptor::detail

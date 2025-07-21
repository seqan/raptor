// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::call_parallel_on_bins.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <algorithm>  // for clamp
#include <bit>        // for bit_ceil
#include <cstddef>    // for size_t
#include <cstdint>    // for uint8_t
#include <functional> // for invoke
#include <ranges>     // for iota_view, __fn, iota, size, views
#include <string>     // for basic_string, string
#include <vector>     // for vector

#include <hibf/contrib/std/zip_view.hpp> // for zip_view, operator-, operator+, operator==, zip, zip_fn
#include <hibf/misc/divide_and_ceil.hpp> // for divide_and_ceil

namespace raptor
{

template <typename algorithm_t>
void call_parallel_on_bins(algorithm_t && worker,
                           std::vector<std::vector<std::string>> const & bin_paths,
                           uint8_t const threads)
{
    size_t const number_of_bins = bin_paths.size();
    // clang-format off
    size_t const chunk_size = std::clamp<size_t>(
        std::bit_ceil(seqan::hibf::divide_and_ceil(number_of_bins, threads)),
        8u,
        64u);
    auto chunked_view = seqan::stl::views::zip(bin_paths, std::views::iota(0u, number_of_bins))
                      | seqan::stl::views::chunk(chunk_size);
    // clang-format on
    size_t const number_of_chunks = std::ranges::size(chunked_view);

#pragma omp parallel for schedule(dynamic) num_threads(threads)
    for (size_t i = 0; i < number_of_chunks; ++i)
    {
        std::invoke(worker, chunked_view[i]);
    }
}

} // namespace raptor

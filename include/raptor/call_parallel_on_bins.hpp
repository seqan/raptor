// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::call_parallel_on_bins.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/algorithm/detail/execution_handler_parallel.hpp>

#include <hibf/contrib/std/chunk_view.hpp>
#include <hibf/contrib/std/zip_view.hpp>

namespace raptor
{

template <typename algorithm_t>
void call_parallel_on_bins(algorithm_t && worker,
                           std::vector<std::vector<std::string>> const & bin_paths,
                           uint8_t const threads)
{
    // GCOVR_EXCL_START
    size_t const chunk_size = std::clamp<size_t>(std::bit_ceil(bin_paths.size() / threads), 8u, 64u);
    // GCOVR_EXCL_STOP
    auto chunked_view = seqan::stl::views::zip(bin_paths, std::views::iota(0u)) | seqan::stl::views::chunk(chunk_size);
    seqan3::detail::execution_handler_parallel executioner{threads};
    executioner.bulk_execute(std::move(worker), std::move(chunked_view), []() {});
}

} // namespace raptor

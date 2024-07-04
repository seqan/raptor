// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::do_parallel.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <functional>
#include <omp.h>
#include <vector>

#include <hibf/misc/divide_and_ceil.hpp>

namespace raptor
{

template <typename algorithm_t>
void do_parallel(algorithm_t && worker, size_t const num_records, size_t const threads)
{
    size_t const chunk_size = seqan::hibf::divide_and_ceil(num_records, threads * threads);
    size_t const number_of_chunks = seqan::hibf::divide_and_ceil(num_records, chunk_size);

#pragma omp parallel for schedule(dynamic) num_threads(threads)
    for (size_t i = 0; i < number_of_chunks; ++i)
    {
        size_t const start = chunk_size * i;
        size_t const extent = i == (number_of_chunks - 1) ? num_records - i * chunk_size : chunk_size;
        std::invoke(worker, start, extent);
    }
}

} // namespace raptor

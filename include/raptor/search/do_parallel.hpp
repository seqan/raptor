// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::do_parallel.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <functional>
#include <omp.h>
#include <vector>

namespace raptor
{

template <typename algorithm_t>
void do_parallel(algorithm_t && worker, size_t const num_records, size_t const threads)
{
    size_t const records_per_thread = num_records / threads;

#pragma omp parallel for schedule(static) num_threads(threads)
    for (size_t i = 0; i < threads; ++i)
    {
        size_t const start = records_per_thread * i;
        size_t const extent = i == (threads - 1) ? num_records - i * records_per_thread : records_per_thread;
        std::invoke(worker, start, extent);
    }
}

} // namespace raptor

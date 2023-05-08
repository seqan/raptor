// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides raptor::do_parallel.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <chrono>
#include <future>
#include <vector>

namespace raptor
{

template <typename algorithm_t>
void do_parallel(algorithm_t && worker, size_t const num_records, size_t const threads)
{
    std::vector<decltype(std::async(std::launch::async, worker, size_t{}, size_t{}))> tasks;
    size_t const records_per_thread = num_records / threads;

    for (size_t i = 0; i < threads; ++i)
    {
        size_t const start = records_per_thread * i;
        size_t const end = i == (threads - 1) ? num_records : records_per_thread * (i + 1);
        tasks.emplace_back(std::async(std::launch::async, worker, start, end));
    }

    for (auto && task : tasks)
        task.get();
}

} // namespace raptor

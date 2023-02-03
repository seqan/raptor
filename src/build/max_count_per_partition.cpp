// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <algorithm>
#include <fstream>

#include <raptor/build/max_count_per_partition.hpp>
#include <raptor/call_parallel_on_bins.hpp>
#include <raptor/file_reader.hpp>

namespace raptor
{

namespace detail
{

template <file_types file_type>
std::vector<size_t> max_count_per_partition(partition_config const & cfg, build_arguments const & arguments)
{
    std::vector<size_t> kmers_per_partition(cfg.partitions);
    std::mutex callback_mutex{};
    file_reader<file_type> const reader{arguments.shape, arguments.window_size};

    auto callback = [&callback_mutex, &kmers_per_partition, &cfg](std::vector<size_t> const & kmer_counts)
    {
        std::lock_guard<std::mutex> guard{callback_mutex};
        for (size_t i = 0; i < cfg.partitions; ++i)
            kmers_per_partition[i] = std::max<size_t>(kmers_per_partition[i], kmer_counts[i]);
    };

    auto worker = [&callback, &reader, &cfg](auto && zipped_view, auto &&)
    {
        std::vector<size_t> max_kmer_counts(cfg.partitions);
        std::vector<size_t> kmer_counts(cfg.partitions);
        for (auto && [file_names, bin_number] : zipped_view)
        {
            reader.for_each_hash(file_names,
                                 [&](auto && hash)
                                 {
                                     ++kmer_counts[cfg.hash_partition(hash)];
                                 });
            for (size_t i = 0; i < cfg.partitions; ++i)
                max_kmer_counts[i] = std::max<size_t>(max_kmer_counts[i], kmer_counts[i]);
            std::ranges::fill(kmer_counts, 0u);
        }
        callback(max_kmer_counts);
    };

    call_parallel_on_bins(worker, arguments.bin_path, arguments.threads);

    return kmers_per_partition;
}

} // namespace detail

std::vector<size_t> max_count_per_partition(partition_config const & cfg, build_arguments const & arguments)
{
    arguments.bin_size_timer.start();
    std::vector<size_t> result = arguments.input_is_minimiser
                                   ? detail::max_count_per_partition<file_types::minimiser>(cfg, arguments)
                                   : detail::max_count_per_partition<file_types::sequence>(cfg, arguments);
    arguments.bin_size_timer.stop();

    return result;
}

} // namespace raptor

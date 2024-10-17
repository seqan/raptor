// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::max_count_per_partition.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <algorithm>
#include <fstream>

#include <hibf/sketch/hyperloglog.hpp>

#include <raptor/build/max_count_per_partition.hpp>
#include <raptor/call_parallel_on_bins.hpp>
#include <raptor/file_reader.hpp>

namespace raptor
{

namespace detail
{

template <file_types file_type>
std::vector<size_t> max_count_per_partition(partition_config const & cfg,
                                            std::vector<std::vector<std::string>> const & bin_path,
                                            uint8_t const threads,
                                            seqan3::shape const & shape,
                                            uint32_t const window_size)
{
    std::vector<size_t> kmers_per_partition(cfg.partitions);
    std::mutex callback_mutex{};
    file_reader<file_type> const reader{shape, window_size};

    auto callback = [&callback_mutex, &kmers_per_partition, &cfg](std::vector<size_t> const & kmer_counts)
    {
        std::lock_guard<std::mutex> guard{callback_mutex};
        for (size_t i = 0; i < cfg.partitions; ++i)
            kmers_per_partition[i] = std::max<size_t>(kmers_per_partition[i], kmer_counts[i]);
    };

    auto worker = [&callback, &reader, &cfg](auto && zipped_view)
    {
        std::vector<size_t> max_kmer_counts(cfg.partitions);
        std::vector<seqan::hibf::sketch::hyperloglog> kmer_counts(cfg.partitions, 15u);

        for (auto && [file_names, bin_number] : zipped_view)
        {
            reader.for_each_hash(file_names,
                                 [&](auto && hash)
                                 {
                                     kmer_counts[cfg.hash_partition(hash)].add(hash);
                                 });
            for (size_t i = 0; i < cfg.partitions; ++i)
            {
                max_kmer_counts[i] = std::max<size_t>(max_kmer_counts[i], kmer_counts[i].estimate());
                kmer_counts[i].reset();
            }
        }
        callback(max_kmer_counts);
    };

    call_parallel_on_bins(worker, bin_path, threads);

    return kmers_per_partition;
}

} // namespace detail

std::vector<size_t> max_count_per_partition(partition_config const & cfg, build_arguments const & arguments)
{
    arguments.bin_size_timer.start();
    // GCOVR_EXCL_START
    std::vector<size_t> result = arguments.input_is_minimiser
                                   ? detail::max_count_per_partition<file_types::minimiser>(cfg,
                                                                                            arguments.bin_path,
                                                                                            arguments.threads,
                                                                                            arguments.shape,
                                                                                            arguments.window_size)
                                   : detail::max_count_per_partition<file_types::sequence>(cfg,
                                                                                           arguments.bin_path,
                                                                                           arguments.threads,
                                                                                           arguments.shape,
                                                                                           arguments.window_size);
    // GCOVR_EXCL_STOP
    arguments.bin_size_timer.stop();

    return result;
}

} // namespace raptor

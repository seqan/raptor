// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <algorithm>
#include <fstream>

#include <raptor/build/max_count_per_partition.hpp>

namespace raptor
{

std::vector<size_t> max_count_per_partition(partition_config const & cfg,
                                            std::vector<std::vector<std::string>> const & bin_path)
{
    std::vector<size_t> kmers_per_partition(cfg.partitions);

    auto callback = [&kmers_per_partition, &cfg](std::vector<size_t> const & kmer_counts)
    {
        for (size_t i = 0; i < cfg.partitions; ++i)
            kmers_per_partition[i] = std::max<size_t>(kmers_per_partition[i], kmer_counts[i]);
    };

    std::vector<size_t> kmer_counts(cfg.partitions);

    for (size_t part = 0; part < cfg.partitions; ++part)
    {
        uint64_t hash;
        for (auto && file_names : bin_path)
        {
            for (auto && file_name : file_names)
            {
                std::ifstream infile{file_name, std::ios::binary};

                while (infile.read(reinterpret_cast<char *>(&hash), sizeof(hash)))
                    ++kmer_counts[(hash & cfg.mask) / cfg.suffixes_per_part];
            }
            callback(kmer_counts);
            std::ranges::fill(kmer_counts, 0u);
        }
    }

    return kmers_per_partition;
}

} // namespace raptor

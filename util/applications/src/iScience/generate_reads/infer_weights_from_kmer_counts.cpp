// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <hibf/contrib/robin_hood.hpp>

#include <raptor/file_reader.hpp>

#include "cmd_arguments.hpp"

namespace raptor::util::generate_reads
{

void infer_weights_from_kmer_counts(cmd_arguments & arguments)
{
    raptor::file_reader<raptor::file_types::sequence> const reader{seqan3::ungapped{32u}, 32u};

#pragma omp parallel for schedule(dynamic) num_threads(arguments.threads)
    for (size_t i = 0; i < arguments.number_of_bins; ++i)
    {
        // could also be defined outside the loop and then be used with the pragma directtive `private(kmer_set)`
        // and then be cleared at the end of the loop
        // However, this would increase the maximum memory usage, as the set would stay at the maximum size of all
        // bins that are processed by this thread.
        robin_hood::unordered_flat_set<uint64_t> kmer_set{};
        reader.hash_into(arguments.bin_path[i], std::inserter(kmer_set, kmer_set.begin()));
        arguments.number_of_reads_per_bin[i] = kmer_set.size();
    }
}

} // namespace raptor::util::generate_reads

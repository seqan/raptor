// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <vector>

namespace raptor::util::generate_reads
{

enum class weight : uint8_t
{
    from_file_sizes,
    from_hll_sketches,
    from_kmer_counts,
    from_weight_column,
    from_uniform_distribution
};

struct cmd_arguments
{
    // User input
    weight weight_mode{weight::from_hll_sketches};
    std::filesystem::path bin_file_path{};
    std::filesystem::path output_filename{};
    uint8_t errors{2u};
    uint32_t read_length{100u};
    uint64_t number_of_reads{1ULL << 20};
    uint64_t threads{1u};
    bool print_weights{false};

    // Internal state
    size_t number_of_bins{};
    std::vector<std::filesystem::path> bin_path{};
    std::vector<size_t> number_of_reads_per_bin{};
};

} // namespace raptor::util::generate_reads

// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::build_arguments.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <cstdint>    // for uint64_t, uint8_t, uint32_t
#include <filesystem> // for path
#include <string>     // for basic_string, string
#include <vector>     // for vector

#include <seqan3/search/kmer_index/shape.hpp> // for shape, ungapped

#include <hibf/misc/timer.hpp> // for concurrent_timer

namespace raptor
{

struct build_arguments
{
    // Related to k-mers
    uint8_t kmer_size{20u};
    uint32_t window_size{kmer_size};
    std::string shape_string{};
    seqan3::shape shape{seqan3::ungapped{kmer_size}};

    // Related to IBF
    std::filesystem::path out_path{};
    uint64_t bins{64};
    mutable uint64_t bits{4096}; // Allow to change bits for each partition
    uint64_t hash{2};
    uint8_t parts{1u};
    double fpr{0.05};

    // General arguments
    std::vector<std::vector<std::string>> bin_path{};
    std::filesystem::path bin_file{};
    uint8_t threads{1u};
    bool is_hibf{false};
    bool input_is_minimiser{false};
    bool quiet{false};
    std::filesystem::path timing_out{};

    // Timers do not copy the stored duration upon copy construction/assignment
    mutable seqan::hibf::concurrent_timer wall_clock_timer{};
    mutable seqan::hibf::concurrent_timer bin_size_timer{};
    mutable seqan::hibf::concurrent_timer index_allocation_timer{};
    mutable seqan::hibf::concurrent_timer user_bin_io_timer{};
    mutable seqan::hibf::concurrent_timer merge_kmers_timer{};
    mutable seqan::hibf::concurrent_timer fill_ibf_timer{};
    mutable seqan::hibf::concurrent_timer store_index_timer{};

    void print_timings() const;
    void write_timings_to_file() const;
};

} // namespace raptor

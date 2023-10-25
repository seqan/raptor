// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides raptor::build_arguments.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <filesystem>
#include <vector>

#include <seqan3/search/kmer_index/shape.hpp>

#include <hibf/misc/timer.hpp>

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

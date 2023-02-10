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

#include <raptor/argument_parsing/memory_usage.hpp>
#include <raptor/argument_parsing/timer.hpp>
#include <raptor/strong_types.hpp>

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
    std::filesystem::path out_path{"./"};
    uint64_t bins{64};
    mutable uint64_t bits{4096}; // Allow to change bits for each partition
    uint64_t hash{2};
    uint8_t parts{1u};
    double fpr{0.05};
    bool compressed{false};

    // General arguments
    std::vector<std::vector<std::string>> bin_path{};
    std::filesystem::path bin_file{};
    uint8_t threads{1u};
    bool is_hibf{false};
    bool input_is_minimiser{false};
    bool verbose{false};

    // Timers do not copy the stored duration upon copy construction/assignment
    mutable timer<concurrent::yes> wall_clock_timer{};
    mutable timer<concurrent::yes> bin_size_timer{};
    mutable timer<concurrent::yes> index_allocation_timer{};
    mutable timer<concurrent::yes> user_bin_io_timer{};
    mutable timer<concurrent::yes> merge_kmers_timer{};
    mutable timer<concurrent::yes> fill_ibf_timer{};
    mutable timer<concurrent::yes> store_index_timer{};

    // GCOVR_EXCL_START
    void print_timings() const
    {
        if (!verbose)
            return;
        std::cerr << std::fixed << std::setprecision(2) << "============= Timings =============\n";
        std::cerr << "Wall clock time [s]: " << wall_clock_timer.in_seconds() << '\n';
        std::cerr << "Peak memory usage " << formatted_peak_ram() << '\n';
        if (!is_hibf)
            std::cerr << "Determine IBF size [s]: " << bin_size_timer.in_seconds() << '\n';
        std::cerr << "Index allocation [s]: " << index_allocation_timer.in_seconds() << '\n';
        std::cerr << "User bin I/O avg per thread [s]: " << user_bin_io_timer.in_seconds() / threads << '\n';
        std::cerr << "User bin I/O sum [s]: " << user_bin_io_timer.in_seconds() << '\n';
        if (is_hibf)
        {
            std::cerr << "Merge kmer sets avg per thread [s]: " << merge_kmers_timer.in_seconds() / threads << '\n';
            std::cerr << "Merge kmer sets sum [s]: " << merge_kmers_timer.in_seconds() << '\n';
        }
        std::cerr << "Fill IBF avg per thread [s]: " << fill_ibf_timer.in_seconds() / threads << '\n';
        std::cerr << "Fill IBF sum [s]: " << fill_ibf_timer.in_seconds() << '\n';
        std::cerr << "Store index [s]: " << store_index_timer.in_seconds() << '\n';
    }
    // GCOVR_EXCL_STOP
};

} // namespace raptor

// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides raptor::search_arguments.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <filesystem>
#include <vector>

#include <seqan3/search/kmer_index/shape.hpp>

#include <raptor/argument_parsing/formatted_index_size.hpp>
#include <raptor/argument_parsing/memory_usage.hpp>
#include <raptor/threshold/threshold_parameters.hpp>

#include <hibf/misc/timer.hpp>

namespace raptor
{

struct search_arguments
{
    // Related to k-mers
    uint32_t window_size{20u};
    seqan3::shape shape{seqan3::ungapped{20u}};
    uint8_t shape_size{shape.size()};
    uint8_t shape_weight{shape.count()};
    uint8_t threads{1u};
    uint8_t parts{1u};

    // Related to thresholding
    double tau{0.9999};
    double threshold{std::numeric_limits<double>::quiet_NaN()};
    double p_max{0.15};
    double fpr{0.05};
    uint64_t query_length{};
    uint8_t errors{0};

    // Related to IBF
    std::filesystem::path index_file{};

    // General arguments
    std::vector<std::vector<std::string>> bin_path{};
    std::filesystem::path query_file{};
    std::filesystem::path out_file{"search.out"};
    bool write_time{false};
    bool is_hibf{false};
    bool cache_thresholds{false};
    bool quiet{false};
    std::filesystem::path timing_out{};

    // Timers do not copy the stored duration upon copy construction/assignment
    mutable seqan::hibf::concurrent_timer wall_clock_timer{};
    mutable seqan::hibf::concurrent_timer query_length_timer{};
    mutable seqan::hibf::concurrent_timer query_file_io_timer{};
    mutable seqan::hibf::concurrent_timer load_index_timer{};
    mutable seqan::hibf::concurrent_timer compute_minimiser_timer{};
    mutable seqan::hibf::concurrent_timer query_ibf_timer{};
    mutable seqan::hibf::concurrent_timer generate_results_timer{};
    mutable seqan::hibf::concurrent_timer complete_search_timer{};
    mutable seqan::hibf::concurrent_timer parallel_search_timer{};

    void print_timings() const;
    void write_timings_to_file() const;

    raptor::threshold::threshold_parameters make_threshold_parameters() const noexcept
    {
        return {.window_size = window_size,
                .shape = shape,
                .query_length = query_length,
                .errors = errors,
                .percentage = threshold,
                .p_max = p_max,
                .tau = tau,
                .cache_thresholds = cache_thresholds,
                .output_directory = index_file.parent_path()};
    }
};

} // namespace raptor

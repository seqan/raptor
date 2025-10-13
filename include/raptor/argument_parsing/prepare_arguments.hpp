// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::prepare_arguments.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <cstdint>    // for uint8_t, uint32_t
#include <filesystem> // for path
#include <iomanip>    // for operator<<, setprecision
#include <iostream>   // for basic_ostream, operator<<, cerr, fixed
#include <string>     // for basic_string, char_traits, string
#include <vector>     // for vector

#include <seqan3/search/kmer_index/shape.hpp> // for shape, ungapped

#include <hibf/misc/timer.hpp> // for concurrent_timer

#include <raptor/argument_parsing/memory_usage.hpp> // for formatted_peak_ram

namespace raptor
{

struct prepare_arguments
{
    uint8_t kmer_size{20u};
    uint32_t window_size{kmer_size};
    std::string shape_string{};
    seqan3::shape shape{seqan3::ungapped{kmer_size}};
    bool use_filesize_dependent_cutoff{false};
    uint8_t kmer_count_cutoff{1u};

    std::filesystem::path out_dir{"./"};

    std::vector<std::vector<std::string>> bin_path{};
    std::filesystem::path bin_file{};
    uint8_t threads{1u};
    bool quiet{false};

    // Timers do not copy the stored duration upon copy construction/assignment
    mutable seqan::hibf::concurrent_timer wall_clock_timer{};
    mutable seqan::hibf::concurrent_timer compute_minimiser_timer{};
    mutable seqan::hibf::concurrent_timer write_minimiser_timer{};
    mutable seqan::hibf::concurrent_timer write_header_timer{};

    void print_timings() const
    {
        if (quiet)
            return;
        std::cerr << std::fixed << std::setprecision(2) << "============= Timings =============\n";
        std::cerr << "Wall clock time [s]: " << wall_clock_timer.in_seconds() << '\n';
        std::cerr << "Peak memory usage " << formatted_peak_ram() << '\n';
        std::cerr << "Compute minimiser [s]: " << compute_minimiser_timer.in_seconds() / threads << '\n';
        std::cerr << "Write minimiser files [s]: " << write_minimiser_timer.in_seconds() / threads << '\n';
        std::cerr << "Write header files [s]: " << write_header_timer.in_seconds() / threads << '\n';
    }
};

} // namespace raptor

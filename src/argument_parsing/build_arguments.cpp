// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Implements raptor::build_arguments.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <fstream>

#include <raptor/argument_parsing/build_arguments.hpp>
#include <raptor/argument_parsing/formatted_index_size.hpp>
#include <raptor/argument_parsing/memory_usage.hpp>

namespace raptor
{

void build_arguments::print_timings() const
{
    std::cerr << std::fixed << std::setprecision(2) << "============= Timings =============\n";
    std::cerr << "Wall clock time [s]: " << wall_clock_timer.in_seconds() << '\n';
    std::cerr << "Peak memory usage " << formatted_peak_ram() << '\n';
    std::cerr << "Index size " << formatted_index_size(out_path, parts) << '\n';

    if (is_hibf)
        std::cerr << "Determine IBF size [s]: NA\n";
    else
        std::cerr << "Determine IBF size [s]: " << bin_size_timer.in_seconds() << '\n';

    std::cerr << "Index allocation [s]: " << index_allocation_timer.in_seconds() << '\n';
    std::cerr << "User bin I/O avg per thread [s]: " << user_bin_io_timer.in_seconds() / threads << '\n';
    std::cerr << "User bin I/O sum [s]: " << user_bin_io_timer.in_seconds() << '\n';

    if (is_hibf)
    {
        std::cerr << "Merge kmer sets avg per thread [s]: " << merge_kmers_timer.in_seconds() / threads << '\n';
        std::cerr << "Merge kmer sets sum [s]: " << merge_kmers_timer.in_seconds() << '\n';
    }
    else
    {
        std::cerr << "Merge kmer sets avg per thread [s]: NA\n";
        std::cerr << "Merge kmer sets sum [s]: NA\n";
    }

    std::cerr << "Fill IBF avg per thread [s]: " << fill_ibf_timer.in_seconds() / threads << '\n';
    std::cerr << "Fill IBF sum [s]: " << fill_ibf_timer.in_seconds() << '\n';
    std::cerr << "Store index [s]: " << store_index_timer.in_seconds() << '\n';
}

void build_arguments::write_timings_to_file() const
{
    std::ofstream output_stream{timing_out};
    output_stream << std::fixed << std::setprecision(2);
    output_stream << "wall_clock_time_in_seconds\t"
                  << "peak_memory_usage_in_kibibytes\t"
                  << "index_size_in_kibibytes\t"
                  << "determine_ibf_size_in_seconds\t"
                  << "index_allocation_in_seconds\t"
                  << "user_bin_io_avg_per_thread_in_seconds\t"
                  << "user_bin_io_sum_in_seconds\t"
                  << "merge_kmer_sets_avg_per_thread_in_seconds\t"
                  << "merge_kmer_sets_sum_in_seconds\t"
                  << "fill_ibf_avg_per_thread_in_seconds\t"
                  << "fill_ibf_sum_in_seconds\t"
                  << "store_index_in_seconds\n";

    output_stream << wall_clock_timer.in_seconds() << '\t';

    if (long const peak_ram_KiB = peak_ram_in_KiB(); peak_ram_KiB != -1L)
        output_stream << peak_ram_KiB << '\t';
    else
        output_stream << "NA\t"; // GCOVR_EXCL_LINE

    output_stream << index_size_in_KiB(out_path, parts) << '\t';

    if (is_hibf)
        output_stream << "NA\t";
    else
        output_stream << bin_size_timer.in_seconds() << '\t';

    output_stream << index_allocation_timer.in_seconds() << '\t';
    output_stream << user_bin_io_timer.in_seconds() / threads << '\t';
    output_stream << user_bin_io_timer.in_seconds() << '\t';

    if (is_hibf)
    {
        output_stream << merge_kmers_timer.in_seconds() / threads << '\t';
        output_stream << merge_kmers_timer.in_seconds() << '\t';
    }
    else
    {
        output_stream << "NA\t";
        output_stream << "NA\t";
    }

    output_stream << fill_ibf_timer.in_seconds() / threads << '\t';
    output_stream << fill_ibf_timer.in_seconds() << '\t';
    output_stream << store_index_timer.in_seconds() << '\n';
}

} // namespace raptor

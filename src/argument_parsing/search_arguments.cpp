// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Implements raptor::search_arguments.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <fstream>

#include <raptor/argument_parsing/formatted_index_size.hpp>
#include <raptor/argument_parsing/memory_usage.hpp>
#include <raptor/argument_parsing/search_arguments.hpp>

namespace raptor
{

void search_arguments::print_timings() const
{
    std::cerr << std::fixed << std::setprecision(2) << "============= Timings =============\n";
    std::cerr << "Wall clock time [s]: " << wall_clock_timer.in_seconds() << '\n';
    std::cerr << "Peak memory usage " << formatted_peak_ram() << '\n';
    std::cerr << "Index size " << formatted_index_size(index_file, parts) << '\n';
    std::cerr << "Determine query length [s]: " << query_length_timer.in_seconds() << '\n';
    std::cerr << "Query file I/O [s]: " << query_file_io_timer.in_seconds() << '\n';
    std::cerr << "Load index [s]: " << load_index_timer.in_seconds() << '\n';
    std::cerr << "Compute minimiser avg per thread [s]: " << compute_minimiser_timer.in_seconds() / threads << '\n';
    std::cerr << "Compute minimiser sum [s]: " << compute_minimiser_timer.in_seconds() << '\n';
    std::cerr << "Query IBF avg per thread [s]: " << query_ibf_timer.in_seconds() / threads << '\n';
    std::cerr << "Query IBF sum [s]: " << query_ibf_timer.in_seconds() << '\n';
    std::cerr << "Generate results avg per thread [s]: " << generate_results_timer.in_seconds() / threads << '\n';
    std::cerr << "Generate results sum [s]: " << generate_results_timer.in_seconds() << '\n';
}

void search_arguments::write_timings_to_file() const
{
    std::ofstream output_stream{timing_out, std::ios_base::app};
    output_stream << std::fixed << std::setprecision(2);
    output_stream << "wall_clock_time_in_seconds\t"
                  << "peak_memory_usage_in_KiB\t"
                  << "index_size_in_KiB\t"
                  << "determine_query_length_in_seconds\t"
                  << "query_file_io_in_seconds\t"
                  << "load_index_in_seconds\t"
                  << "compute_minimiser_avg_per_thread_in_seconds\t"
                  << "compute_minimiser_sum_in_seconds\t"
                  << "query_ibf_avg_per_thread_in_seconds\t"
                  << "query_ibf_sum_in_seconds\t"
                  << "generate_results_avg_per_thread_in_seconds\t"
                  << "generate_results_sum_in_seconds\n";

    output_stream << wall_clock_timer.in_seconds() << '\t';

    if (long const peak_ram_KiB = peak_ram_in_KiB(); peak_ram_KiB != -1L)
        output_stream << peak_ram_KiB << '\t';
    else
        output_stream << "NA\t"; // GCOVR_EXCL_LINE

    output_stream << index_size_in_KiB(index_file, parts) << '\t';
    output_stream << query_length_timer.in_seconds() << '\t';
    output_stream << query_file_io_timer.in_seconds() << '\t';
    output_stream << load_index_timer.in_seconds() << '\t';
    output_stream << compute_minimiser_timer.in_seconds() / threads << '\t';
    output_stream << compute_minimiser_timer.in_seconds() << '\t';
    output_stream << query_ibf_timer.in_seconds() / threads << '\t';
    output_stream << query_ibf_timer.in_seconds() << '\t';
    output_stream << generate_results_timer.in_seconds() / threads << '\t';
    output_stream << generate_results_timer.in_seconds() << '\n';
}

} // namespace raptor

// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::search_arguments.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <fstream>

#include <raptor/argument_parsing/cpu_time.hpp>
#include <raptor/argument_parsing/formatted_index_size.hpp>
#include <raptor/argument_parsing/memory_usage.hpp>
#include <raptor/argument_parsing/search_arguments.hpp>

namespace raptor
{

void search_arguments::print_timings() const
{
    std::cerr << std::fixed << std::setprecision(2) << "============= Timings =============\n";
    std::cerr << "Peak memory usage " << formatted_peak_ram() << '\n';
    std::cerr << "Index size " << formatted_index_size(index_file, parts) << '\n';
    std::cerr << "Configured threads: " << static_cast<size_t>(threads) << '\n';
    std::cerr << "Wall clock time [s]: " << wall_clock_timer.in_seconds() << '\n';

    double cpu_usage_search{-1.0};
    if (cpu_time_t cpu_time = get_cpu_time(); cpu_time.is_valid())
    {
        std::cerr << "├── User time [s]: " << cpu_time.user_time_in_seconds << '\n';
        std::cerr << "├── System time [s]: " << cpu_time.system_time_in_seconds << '\n';
        std::cerr << "├── CPU usage [%]: " << cpu_time.cpu_usage_in_percent(wall_clock_timer.in_seconds(), threads)
                  << '\n';

        // Substract user/system times that do not belong to the parallel search.
        cpu_time -= (wall_clock_timer.in_seconds() - parallel_search_timer.in_seconds());
        if (cpu_time.is_valid())
            cpu_usage_search = cpu_time.cpu_usage_in_percent(parallel_search_timer.in_seconds(), threads);
    }
    else
    {
        std::cerr << "├── User time [s]: Not available\n";   // GCOVR_EXCL_LINE
        std::cerr << "├── System time [s]: Not available\n"; // GCOVR_EXCL_LINE
        std::cerr << "├── CPU usage [%]: Not available\n";   // GCOVR_EXCL_LINE
    }

    std::cerr << "├── Determine query length [s]: " << query_length_timer.in_seconds() << '\n';
    std::cerr << "└── Complete search [s]: " << complete_search_timer.in_seconds() << '\n';
    std::cerr << "    ├── Query file I/O [s]: " << query_file_io_timer.in_seconds() << '\n';
    std::cerr << "    ├── Load index [s]: " << load_index_timer.in_seconds() << '\n';
    std::cerr << "    └── Parallel search [s]: " << parallel_search_timer.in_seconds() << '\n';

    if (cpu_usage_search > 0.0)
        std::cerr << "        ├── CPU usage [%]: " << cpu_usage_search << '\n';
    else
        std::cerr << "        ├── CPU usage [%]: Not available\n"; // GCOVR_EXCL_LINE

    std::cerr << "        ├── Compute minimiser\n";
    std::cerr << "        │   ├── Max [s]: " << compute_minimiser_timer.max_in_seconds() << '\n';
    std::cerr << "        │   └── Avg [s]: " << compute_minimiser_timer.avg_in_seconds() << '\n';
    std::cerr << "        ├── Query IBF\n";
    std::cerr << "        │   ├── Max [s]: " << query_ibf_timer.max_in_seconds() << '\n';
    std::cerr << "        │   └── Avg [s]: " << query_ibf_timer.avg_in_seconds() << '\n';
    std::cerr << "        └── Generate results\n";
    std::cerr << "            ├── Max [s]: " << generate_results_timer.max_in_seconds() << '\n';
    std::cerr << "            └── Avg [s]: " << generate_results_timer.avg_in_seconds() << '\n';
}

void search_arguments::write_timings_to_file() const
{
    std::ofstream output_stream{timing_out, std::ios_base::app};
    output_stream << std::fixed << std::setprecision(2);
    output_stream << "peak_memory_usage_in_kibibytes\t"
                  << "index_size_in_kibibytes\t"
                  << "configured_threads\t"
                  << "wall_clock_time_in_seconds\t"
                  << "user_time_in_seconds\t"
                  << "system_time_in_seconds\t"
                  << "cpu_usage_in_percent\t"
                  << "determine_query_length_in_seconds\t"
                  << "complete_search_in_seconds\t"
                  << "query_file_io_in_seconds\t"
                  << "load_index_in_seconds\t"
                  << "parallel_search_in_seconds\t"
                  << "cpu_usage_parallel_search_in_percent\t"
                  << "compute_minimiser_max_in_seconds\t"
                  << "compute_minimiser_avg_in_seconds\t"
                  << "query_ibf_max_in_seconds\t"
                  << "query_ibf_avg_in_seconds\t"
                  << "generate_results_max_in_seconds\t"
                  << "generate_results_avg_in_seconds\n";

    if (long const peak_ram_KiB = peak_ram_in_KiB(); peak_ram_KiB != -1L)
        output_stream << peak_ram_KiB << '\t';
    else
        output_stream << "NA\t"; // GCOVR_EXCL_LINE

    output_stream << index_size_in_KiB(index_file, parts) << '\t';
    output_stream << static_cast<size_t>(threads) << '\t';
    output_stream << wall_clock_timer.in_seconds() << '\t';

    double cpu_usage_search{-1.0};
    if (cpu_time_t cpu_time = get_cpu_time(); cpu_time.is_valid())
    {
        output_stream << cpu_time.user_time_in_seconds << '\t';
        output_stream << cpu_time.system_time_in_seconds << '\t';
        output_stream << cpu_time.cpu_usage_in_percent(wall_clock_timer.in_seconds(), threads) << '\t';

        // Substract user/system times that do not belong to the parallel search.
        cpu_time -= (wall_clock_timer.in_seconds() - parallel_search_timer.in_seconds());
        if (cpu_time.is_valid())
            cpu_usage_search = cpu_time.cpu_usage_in_percent(parallel_search_timer.in_seconds(), threads);
    }
    else
    {
        output_stream << "NA\t"; // GCOVR_EXCL_LINE
        output_stream << "NA\t"; // GCOVR_EXCL_LINE
        output_stream << "NA\t"; // GCOVR_EXCL_LINE
    }

    output_stream << query_length_timer.in_seconds() << '\t';
    output_stream << complete_search_timer.in_seconds() << '\t';
    output_stream << query_file_io_timer.in_seconds() << '\t';
    output_stream << load_index_timer.in_seconds() << '\t';
    output_stream << parallel_search_timer.in_seconds() << '\t';

    if (cpu_usage_search > 0.0)
        output_stream << cpu_usage_search << '\t';
    else
        output_stream << "NA\t"; // GCOVR_EXCL_LINE

    output_stream << compute_minimiser_timer.max_in_seconds() << '\t';
    output_stream << compute_minimiser_timer.avg_in_seconds() << '\t';
    output_stream << query_ibf_timer.max_in_seconds() << '\t';
    output_stream << query_ibf_timer.avg_in_seconds() << '\t';
    output_stream << generate_results_timer.max_in_seconds() << '\t';
    output_stream << generate_results_timer.avg_in_seconds() << '\n';
}

} // namespace raptor

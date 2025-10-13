// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::compute_minimiser.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <algorithm>  // for min
#include <cstddef>    // for size_t
#include <cstdint>    // for uint8_t, uint64_t, uint16_t
#include <filesystem> // for path, exists, remove
#include <fstream>    // for basic_ofstream, basic_ostream, operator<<, basic_ios
#include <functional> // for invoke
#include <ranges>     // for iota_view, __fn, operator==, views, iota, size
#include <string>     // for basic_string, char_traits
#include <vector>     // for vector

#include <seqan3/contrib/std/chunk_view.hpp>           // for chunk_view, chunk, chunk_fn
#include <seqan3/contrib/std/detail/adaptor_base.hpp>  // for operator|
#include <seqan3/contrib/std/pair.hpp>                 // for get
#include <seqan3/contrib/std/zip_view.hpp>             // for zip_view, operator-, operator+, operator==, zip
#include <seqan3/search/kmer_index/shape.hpp>          // for shape
#include <seqan3/search/views/kmer_hash.hpp>           // for operator-, operator==
#include <seqan3/search/views/minimiser.hpp>           // for operator!=
#include <seqan3/utility/container/dynamic_bitset.hpp> // for operator==

#include <hibf/contrib/robin_hood.hpp>   // for pair, unordered_map
#include <hibf/misc/divide_and_ceil.hpp> // for divide_and_ceil
#include <hibf/misc/timer.hpp>           // for serial_timer, concurrent_timer

#include <raptor/argument_parsing/prepare_arguments.hpp> // for prepare_arguments
#include <raptor/file_reader.hpp>                        // for file_types, file_reader
#include <raptor/prepare/compute_minimiser.hpp>          // for compute_minimiser
#include <raptor/prepare/cutoff.hpp>                     // for cutoff

namespace raptor
{

std::filesystem::path get_output_path(std::filesystem::path const & output_dir, std::filesystem::path const & file_name)
{
    std::filesystem::path result{output_dir};
    bool const is_compressed = raptor::cutoff::file_is_compressed(file_name);
    result /= is_compressed ? file_name.stem().stem() : file_name.stem();
    result += ".dummy_extension"; // https://github.com/seqan/raptor/issues/355
    return result;
}

void write_list_file(prepare_arguments const & arguments)
{
    std::filesystem::path list_file = arguments.out_dir;
    list_file /= "minimiser.list";
    std::ofstream file{list_file};

    for (auto && file_names : arguments.bin_path)
    {
        std::filesystem::path file_path = get_output_path(arguments.out_dir, file_names[0]);
        file_path.replace_extension("minimiser");
        file << file_path.c_str() << '\n';
    }
}

void compute_minimiser(prepare_arguments const & arguments)
{
    file_reader<file_types::sequence> const reader{arguments.shape, arguments.window_size};
    raptor::cutoff const cutoffs{arguments};

    auto worker = [&](auto && zipped_view)
    {
        seqan::hibf::serial_timer local_compute_minimiser_timer{};
        seqan::hibf::serial_timer local_write_minimiser_timer{};
        seqan::hibf::serial_timer local_write_header_timer{};

        for (auto && [file_names, bin_number] : zipped_view)
        {
            std::filesystem::path const file_name{file_names[0]};
            std::filesystem::path output_path = get_output_path(arguments.out_dir, file_name);

            std::filesystem::path const minimiser_file =
                std::filesystem::path{output_path}.replace_extension("minimiser");
            std::filesystem::path const progress_file =
                std::filesystem::path{output_path}.replace_extension("in_progress");
            std::filesystem::path const header_file = std::filesystem::path{output_path}.replace_extension("header");

            // If we are already done with this file, we can skip it. Otherwise, we create a ".in_progress" file to keep
            // track of whether the minimiser computation was successful.
            bool const already_done = std::filesystem::exists(minimiser_file) && std::filesystem::exists(header_file)
                                   && !std::filesystem::exists(progress_file);

            if (already_done)
                continue;
            else
                std::ofstream outfile{progress_file, std::ios::binary};

            // The hash table stores how often a minimiser appears. It does not matter whether a minimiser appears
            // 50 times or 2000 times, it is stored regardless because the biggest cutoff value is 50. Hence,
            // the hash table stores only values up to 254 to save memory.
            robin_hood::unordered_map<uint64_t, uint8_t> minimiser_table{};
            // The map is (re-)constructed for each file. The alternative is to construct it once for each thread
            // and clear+reuse it for every file that a thread works on. However, this dramatically increases
            // memory consumption because the map will stay as big as needed for the biggest encountered file.

            local_compute_minimiser_timer.start();
            reader.for_each_hash(file_names,
                                 [&](auto && hash)
                                 {
                                     minimiser_table[hash] = std::min<uint8_t>(254u, minimiser_table[hash] + 1);
                                 });
            local_compute_minimiser_timer.stop();

            uint8_t const cutoff = cutoffs.get(file_name);
            uint64_t count{};

            local_write_minimiser_timer.start();
            {
                std::ofstream outfile{minimiser_file, std::ios::binary};
                for (auto && [hash, occurrences] : minimiser_table)
                {
                    if (occurrences >= cutoff)
                    {
                        outfile.write(reinterpret_cast<char const *>(&hash), sizeof(hash));
                        ++count;
                    }
                }
            }
            local_write_minimiser_timer.stop();

            local_write_header_timer.start();
            {
                std::ofstream headerfile{header_file};
                headerfile << arguments.shape.to_string() << '\t' << arguments.window_size << '\t'
                           << static_cast<uint16_t>(cutoff) << '\t' << count << '\n';
            }
            local_write_header_timer.stop();

            std::filesystem::remove(progress_file);
        }

        arguments.compute_minimiser_timer += local_compute_minimiser_timer;
        arguments.write_minimiser_timer += local_write_minimiser_timer;
        arguments.write_header_timer += local_write_header_timer;
    };

    size_t const number_of_bins = arguments.bin_path.size();
    size_t const chunk_size = seqan::hibf::divide_and_ceil(number_of_bins, arguments.threads);
    auto chunked_view = seqan::stl::views::zip(arguments.bin_path, std::views::iota(0u, number_of_bins))
                      | seqan::stl::views::chunk(chunk_size);
    size_t const number_of_chunks = std::ranges::size(chunked_view);

#pragma omp parallel for schedule(dynamic) num_threads(arguments.threads)
    for (size_t i = 0; i < number_of_chunks; ++i)
    {
        std::invoke(worker, chunked_view[i]);
    }

    write_list_file(arguments);
}

} // namespace raptor

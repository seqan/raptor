// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::compute_bin_size and raptor::max_bin_count.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <cassert>    // for assert
#include <cstddef>    // for size_t
#include <cstdint>    // for uint64_t, uint8_t, uint16_t, uint32_t
#include <filesystem> // for path, file_size
#include <fstream>    // for basic_ifstream, basic_istream, operator>>, basic_ios
#include <functional> // for equal_to
#include <iterator>   // for insert_iterator, inserter
#include <mutex>      // for mutex, lock_guard
#include <ranges>     // for operator==, views
#include <string>     // for basic_string, char_traits, string
#include <utility>    // for move
#include <vector>     // for vector

#include <seqan3/search/kmer_index/shape.hpp>          // for shape
#include <seqan3/search/views/kmer_hash.hpp>           // for operator-, operator==
#include <seqan3/search/views/minimiser.hpp>           // for operator!=
#include <seqan3/utility/container/dynamic_bitset.hpp> // for operator==

#include <hibf/build/bin_size_in_bits.hpp> // for bin_size_in_bits
#include <hibf/contrib/robin_hood.hpp>     // for Table, hash, unordered_flat_set
#include <hibf/contrib/std/pair.hpp>       // for get
#include <hibf/misc/timer.hpp>             // for concurrent_timer
#include <hibf/sketch/hyperloglog.hpp>     // for hyperloglog

#include <raptor/argument_parsing/build_arguments.hpp>  // for build_arguments
#include <raptor/argument_parsing/compute_bin_size.hpp> // for compute_bin_size
#include <raptor/call_parallel_on_bins.hpp>             // for call_parallel_on_bins
#include <raptor/file_reader.hpp>                       // for file_types, file_reader

namespace raptor
{

namespace detail
{

size_t kmer_count_from_minimiser_files(std::vector<std::vector<std::string>> const & bin_path, uint8_t const threads)
{
    std::mutex callback_mutex{};
    size_t max_filesize{};
    std::filesystem::path biggest_file{};

    auto callback = [&callback_mutex, &biggest_file, &max_filesize](std::filesystem::path path, size_t const size)
    {
        std::lock_guard<std::mutex> guard{callback_mutex};
        if (size > max_filesize)
        {
            max_filesize = size;
            biggest_file = std::move(path);
        }
    };

    auto worker = [&callback](auto && zipped_view)
    {
        std::filesystem::path minimiser_file{};
        std::filesystem::path biggest_file{};
        size_t max_filesize{};

        for (auto && [file_names, bin_number] : zipped_view)
        {
            for (auto && file_name : file_names)
            {
                minimiser_file = file_name;
                size_t const size = std::filesystem::file_size(minimiser_file);
                if (size > max_filesize)
                {
                    max_filesize = size;
                    biggest_file = minimiser_file;
                }
            }
        }

        callback(std::move(biggest_file), max_filesize);
    };

    call_parallel_on_bins(worker, bin_path, threads);

    std::string shape_string{};
    uint64_t window_size{};
    uint16_t cutoff{};
    size_t max_count{};

    biggest_file.replace_extension("header");
    std::ifstream file_stream{biggest_file};
    file_stream >> shape_string >> window_size >> cutoff >> max_count;

    return max_count;
}

size_t kmer_count_from_sequence_files(std::vector<std::vector<std::string>> const & bin_path,
                                      uint8_t const threads,
                                      seqan3::shape const & shape,
                                      uint32_t const window_size)
{
    size_t max_count{};
    size_t max_bin_id{};
    std::mutex callback_mutex{};
    file_reader<file_types::sequence> const reader{shape, window_size};

    auto callback = [&](size_t const count, size_t const bin_id)
    {
        std::lock_guard<std::mutex> guard{callback_mutex};
        if (count > max_count)
        {
            max_count = count;
            max_bin_id = bin_id;
        }
    };

    auto worker = [&callback, &reader](auto && zipped_view)
    {
        seqan::hibf::sketch::hyperloglog sketch{15u};

        for (auto && [file_names, bin_number] : zipped_view)
        {
            sketch.reset();
            reader.for_each_hash(file_names,
                                 [&sketch](auto && hash)
                                 {
                                     sketch.add(hash);
                                 });
            callback(sketch.estimate(), bin_number);
        }
    };

    // Use sketches to determine biggest bin.
    call_parallel_on_bins(worker, bin_path, threads);

    // Get exact count for biggest bin. Sketch estimate's accuracy depends on sketch_bits (here: 15).
    robin_hood::unordered_flat_set<uint64_t> kmers{};
    auto insert_it = std::inserter(kmers, kmers.end());
    reader.hash_into(bin_path[max_bin_id], insert_it);

    return kmers.size();
}

} // namespace detail

size_t compute_bin_size(build_arguments const & arguments)
{
    arguments.bin_size_timer.start();
    size_t const max_count = arguments.input_is_minimiser
                               ? detail::kmer_count_from_minimiser_files(arguments.bin_path, arguments.threads)
                               : detail::kmer_count_from_sequence_files(arguments.bin_path,
                                                                        arguments.threads,
                                                                        arguments.shape,
                                                                        arguments.window_size);
    arguments.bin_size_timer.stop();

    assert(max_count > 0u);

    return seqan::hibf::build::bin_size_in_bits(
        {.fpr = arguments.fpr, .hash_count = arguments.hash, .elements = max_count});
}

} // namespace raptor

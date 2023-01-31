// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <seqan3/search/views/minimiser_hash.hpp>

#include <raptor/adjust_seed.hpp>
#include <raptor/argument_parsing/compute_bin_size.hpp>
#include <raptor/build/hibf/bin_size_in_bits.hpp>
#include <raptor/call_parallel_on_bins.hpp>
#include <raptor/dna4_traits.hpp>
#include <raptor/file_reader.hpp>

namespace raptor
{

namespace detail
{

inline size_t kmer_count_from_minimiser_files(raptor::build_arguments const & arguments)
{
    std::mutex callback_mutex{};
    size_t max_filesize{};
    std::filesystem::path biggest_file{};

    auto callback = [&callback_mutex, &biggest_file, &max_filesize](std::filesystem::path const path, size_t const size)
    {
        std::lock_guard<std::mutex> guard{callback_mutex};
        if (size > max_filesize)
        {
            max_filesize = size;
            biggest_file = path;
        }
    };

    auto worker = [&callback](auto && zipped_view, auto &&)
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

        callback(biggest_file, max_filesize);
    };

    call_parallel_on_bins(worker, arguments.bin_path, arguments.threads);

    std::string shape_string{};
    uint64_t window_size{};
    uint16_t cutoff{};
    size_t max_count{};

    biggest_file.replace_extension("header");
    std::ifstream file_stream{biggest_file};
    file_stream >> shape_string >> window_size >> cutoff >> max_count;

    return max_count;
}

inline size_t kmer_count_from_sequence_files(raptor::build_arguments const & arguments)
{
    size_t max_count{};
    std::mutex callback_mutex{};
    file_reader<file_types::sequence> const reader{arguments.shape, arguments.window_size};

    auto callback = [&callback_mutex, &max_count](auto && view)
    {
        auto const count = std::ranges::distance(view);
        {
            std::lock_guard<std::mutex> guard{callback_mutex};
            max_count = std::max<size_t>(max_count, count);
        }
    };

    auto worker = [&callback, &reader](auto && zipped_view, auto &&)
    {
        for (auto && [file_names, bin_number] : zipped_view)
            reader.on_hash(file_names, callback);
    };

    call_parallel_on_bins(worker, arguments.bin_path, arguments.threads);

    return max_count;
}

} // namespace detail

size_t compute_bin_size(raptor::build_arguments const & arguments)
{
    size_t const max_count = arguments.input_is_minimiser ? detail::kmer_count_from_minimiser_files(arguments)
                                                          : detail::kmer_count_from_sequence_files(arguments);

    assert(max_count > 0u);

    return hibf::bin_size_in_bits(arguments, max_count);
}

} // namespace raptor

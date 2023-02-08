// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <seqan3/search/views/minimiser_hash.hpp>

#include <raptor/adjust_seed.hpp>
#include <raptor/build/hibf/insert_into_ibf.hpp>
#include <raptor/file_reader.hpp>

namespace raptor::hibf
{

// automatically does naive splitting if number_of_bins > 1
void insert_into_ibf(robin_hood::unordered_flat_set<size_t> & parent_kmers,
                     robin_hood::unordered_flat_set<size_t> const & kmers,
                     size_t const number_of_bins,
                     size_t const bin_index,
                     seqan3::interleaved_bloom_filter<> & ibf,
                     bool is_root,
                     timer & fill_ibf_timer,
                     timer & merge_kmers_timer)
{
    size_t const chunk_size = kmers.size() / number_of_bins + 1;
    size_t chunk_number{};

    timer local_fill_ibf_timer{};
    local_fill_ibf_timer.start();
    for (auto chunk : kmers | seqan3::views::chunk(chunk_size))
    {
        assert(chunk_number < number_of_bins);
        seqan3::bin_index const bin_idx{bin_index + chunk_number};
        ++chunk_number;
        for (size_t const value : chunk)
            ibf.emplace(value, bin_idx);
    }
    local_fill_ibf_timer.stop();
    fill_ibf_timer += local_fill_ibf_timer;

    timer local_merge_kmers_timer{};
    local_merge_kmers_timer.start();
    if (!is_root)
        parent_kmers.insert(kmers.begin(), kmers.end());
    local_merge_kmers_timer.stop();
    merge_kmers_timer += local_merge_kmers_timer;
}

void insert_into_ibf(build_arguments const & arguments,
                     chopper_pack_record const & record,
                     seqan3::interleaved_bloom_filter<> & ibf)
{
    auto const bin_index = seqan3::bin_index{static_cast<size_t>(record.bin_indices.back())};
    std::vector<uint64_t> values;

    timer local_user_bin_io_timer{};
    local_user_bin_io_timer.start();
    if (arguments.input_is_minimiser)
    {
        file_reader<file_types::minimiser> const reader{};
        reader.hash_into(record.filenames, std::back_inserter(values));
    }
    else
    {
        file_reader<file_types::sequence> const reader{arguments.shape, arguments.window_size};
        reader.hash_into(record.filenames, std::back_inserter(values));
    }
    local_user_bin_io_timer.stop();
    arguments.user_bin_io_timer += local_user_bin_io_timer;

    timer local_fill_ibf_timer{};
    local_fill_ibf_timer.start();
    for (auto && value : values)
        ibf.emplace(value, bin_index);
    local_fill_ibf_timer.stop();
    arguments.fill_ibf_timer += local_fill_ibf_timer;
}

} // namespace raptor::hibf

// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

#include <raptor/search/compute_simple_model.hpp>
#include <raptor/search/do_parallel.hpp>
#include <raptor/search/load_index.hpp>
#include <raptor/search/sync_out.hpp>

namespace raptor
{

template <bool compressed>
void run_program_single_socks(search_arguments const & arguments)
{
    constexpr seqan3::data_layout data_layout_mode = compressed ? seqan3::data_layout::compressed :
                                                                 seqan3::data_layout::uncompressed;
    auto index = raptor_index<data_layout_mode>{};

    double index_io_time{0.0};
    double reads_io_time{0.0};
    double compute_time{0.0};

    auto cereal_worker = [&] ()
    {
        load_index(index, arguments, index_io_time);
    };
    auto cereal_handle = std::async(std::launch::async, cereal_worker);

    std::vector<std::vector<seqan3::dna4>> records{};

    std::ifstream fin{arguments.query_file};

    sync_out synced_out{arguments.out_file};

    auto worker = [&] (size_t const start, size_t const end)
    {
        auto & ibf = index.ibf();
        auto counter = ibf.template counting_agent<uint8_t>();
        std::string result_string{};

        auto hash_view = seqan3::views::minimiser_hash(arguments.shape,
                                                       seqan3::window_size{arguments.window_size},
                                                       seqan3::seed{adjust_seed(arguments.shape_weight)});

        for (auto && seq : records | seqan3::views::slice(start, end))
        {
            result_string.clear();
            for (auto const & elem : seq)
                result_string += seqan3::to_char(elem);
            result_string += ": ";

            auto & result = counter.bulk_count(seq | hash_view);

            constexpr int8_t int_to_char_offset{'0'}; // ASCII offset (usually 48), std::to_string is slow
            for (auto const & elem : result)
                result_string += elem + int_to_char_offset;

            result_string += '\n';
            synced_out.write(result_string);
        }
    };

    std::string line{};
    while (!fin.eof())
    {
        records.clear();
        size_t entries{};

        auto start = std::chrono::high_resolution_clock::now();
        while (entries < (1ULL<<20)*10 && std::getline(fin, line))
        {
            auto v = line | seqan3::views::char_to<seqan3::dna4>;
            records.emplace_back(v.begin(), v.end());
            ++entries;
        }
        auto end = std::chrono::high_resolution_clock::now();
        reads_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();

        cereal_handle.wait();

        do_parallel(worker, records.size(), arguments.threads, compute_time);
    }

// LCOV_EXCL_START
    if (arguments.write_time)
    {
        std::filesystem::path file_path{arguments.out_file};
        file_path += ".time";
        std::ofstream file_handle{file_path};
        file_handle << "Index I/O\tReads I/O\tCompute\n";
        file_handle << std::fixed
                    << std::setprecision(2)
                    << index_io_time << '\t'
                    << reads_io_time << '\t'
                    << compute_time;
    }
// LCOV_EXCL_STOP
}

} // namespace raptor

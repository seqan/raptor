// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <robin_hood.h>

#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

#include <raptor/adjust_seed.hpp>
#include <raptor/call_parallel_on_bins.hpp>
#include <raptor/dna4_traits.hpp>
#include <raptor/file_reader.hpp>
#include <raptor/prepare/compute_minimiser.hpp>
#include <raptor/prepare/cutoff.hpp>

namespace raptor
{

void compute_minimiser(prepare_arguments const & arguments)
{
    file_reader<file_types::sequence> const reader{arguments.shape, arguments.window_size};
    raptor::cutoff const cutoffs{arguments};

    auto worker = [&](auto && zipped_view, auto &&)
    {
        timer local_compute_minimiser_timer{};
        timer local_write_minimiser_timer{};
        timer local_write_header_timer{};

        // The hash table stores how often a minimiser appears. It does not matter whether a minimiser appears
        // 50 times or 2000 times, it is stored regardless because the biggest cutoff value is 50. Hence,
        // the hash table stores only values up to 254 to save memory.
        robin_hood::unordered_map<uint64_t, uint8_t> minimiser_table{};

        for (auto && [file_names, bin_number] : zipped_view)
        {
            local_compute_minimiser_timer.start();
            reader.for_each_hash(file_names,
                                 [&](auto && hash)
                                 {
                                     minimiser_table[hash] = std::min<uint8_t>(254u, minimiser_table[hash] + 1);
                                 });
            local_compute_minimiser_timer.stop();

            std::filesystem::path const file_name{file_names[0]};
            bool const is_compressed = raptor::cutoff::file_is_compressed(file_name);
            uint8_t const cutoff = cutoffs.get(file_name);
            uint64_t count{};

            // Store binary file
            std::filesystem::path output_path{arguments.out_dir};
            output_path /= is_compressed ? file_name.stem().stem() : file_name.stem();
            output_path += ".minimiser";

            local_write_minimiser_timer.start();
            {
                std::ofstream outfile{output_path, std::ios::binary};
                for (auto && [hash, occurrences] : minimiser_table)
                {
                    if (occurrences >= cutoff)
                    {
                        outfile.write(reinterpret_cast<const char *>(&hash), sizeof(hash));
                        ++count;
                    }
                }
            }
            local_write_minimiser_timer.stop();

            // Store header file
            output_path.replace_extension("header");

            local_write_header_timer.start();
            {
                std::ofstream headerfile{output_path};
                headerfile << arguments.shape.to_string() << '\t' << arguments.window_size << '\t'
                           << static_cast<uint16_t>(cutoff) << '\t' << count << '\n';
            }
            local_write_header_timer.stop();

            minimiser_table.clear();
        }

        arguments.compute_minimiser_timer += local_compute_minimiser_timer;
        arguments.write_minimiser_timer += local_write_minimiser_timer;
        arguments.write_header_timer += local_write_header_timer;
    };

    call_parallel_on_bins(worker, arguments.bin_path, arguments.threads);
}

} // namespace raptor

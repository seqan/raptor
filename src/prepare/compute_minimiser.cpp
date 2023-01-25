// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <robin_hood.h>

#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

#include <raptor/adjust_seed.hpp>
#include <raptor/call_parallel_on_bins.hpp>
#include <raptor/dna4_traits.hpp>
#include <raptor/prepare/compute_minimiser.hpp>
#include <raptor/prepare/cutoff.hpp>

namespace raptor
{

void compute_minimiser(prepare_arguments const & arguments)
{
    auto minimiser_view = seqan3::views::minimiser_hash(arguments.shape,
                                                        seqan3::window_size{arguments.window_size},
                                                        seqan3::seed{adjust_seed(arguments.shape.count())});

    raptor::cutoff const cutoffs{arguments};

    auto worker = [&](auto && zipped_view, auto &&)
    {
        // The hash table stores how often a minimiser appears. It does not matter whether a minimiser appears
        // 50 times or 2000 times, it is stored regardless because the biggest cutoff value is 50. Hence,
        // the hash table stores only values up to 254 to save memory.
        robin_hood::unordered_map<uint64_t, uint8_t> minimiser_table{};

        for (auto && [file_names, bin_number] : zipped_view)
        {
            for (auto && file_name : file_names)
            {
                seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>> fin{file_name};

                for (auto & [seq] : fin)
                    for (auto && hash : seq | minimiser_view)
                        minimiser_table[hash] = std::min<uint8_t>(254u, minimiser_table[hash] + 1);
            }

            std::filesystem::path const file_name{file_names[0]};
            bool const is_compressed = raptor::cutoff::file_is_compressed(file_name);
            uint16_t cutoff = cutoffs.get(file_name);
            uint64_t count{};

            // Store binary file
            std::filesystem::path output_path{arguments.out_dir};
            output_path /= is_compressed ? file_name.stem().stem() : file_name.stem();
            output_path += ".minimiser";
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

            // Store header file
            output_path.replace_extension("header");
            {
                std::ofstream headerfile{output_path};
                headerfile << arguments.shape.to_string() << '\t' << arguments.window_size << '\t' << cutoff << '\t'
                           << count << '\n';
            }

            minimiser_table.clear();
        }
    };

    call_parallel_on_bins(worker, arguments.bin_path, arguments.threads);
}

} // namespace raptor

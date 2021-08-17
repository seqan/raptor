// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <raptor/build/store_index.hpp>
#include <raptor/shared.hpp>

namespace raptor
{

template <bool compressed>
void build_from_minimiser(build_arguments const & arguments)
{
    seqan3::interleaved_bloom_filter<> ibf{seqan3::bin_count{arguments.bins},
                                           seqan3::bin_size{arguments.bits / arguments.parts},
                                           seqan3::hash_function_count{arguments.hash}};

    auto worker = [&] (auto && zipped_view, auto &&)
        {
            uint64_t read_number;

            for (auto && [file_names, bin_number] : zipped_view)
            {
                for (auto && file_name : file_names)
                {
                    std::ifstream infile{file_name, std::ios::binary};

                    while(infile.read(reinterpret_cast<char*>(&read_number), sizeof(read_number)))
                        ibf.emplace(read_number, seqan3::bin_index{bin_number});
                }
            }
        };

    call_parallel_on_bins(std::move(worker), arguments);

    if constexpr (compressed)
    {
        seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> cibf{std::move(ibf)};
        store_index(arguments.out_path, cibf, arguments);
    }
    else
    {
        store_index(arguments.out_path, ibf, arguments);
    }
}

} // namespace raptor

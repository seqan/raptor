// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <raptor/build/build_from_minimiser.hpp>
#include <raptor/build/call_parallel_on_bins.hpp>
#include <raptor/build/store_index.hpp>

namespace raptor
{

void build_from_minimiser(build_arguments const & arguments)
{
    raptor_index<> index{arguments};

    auto worker = [&] (auto && zipped_view, auto &&)
        {
            uint64_t read_number;
            auto & ibf = index.ibf();

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

    if (arguments.compressed)
    {
        raptor_index<index_structure::ibf_compressed> cindex{std::move(index)};
        store_index(arguments.out_path, cindex, arguments);
    }
    else
    {
        store_index(arguments.out_path, index, arguments);
    }
}

} // namespace raptor

// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <raptor/build/hibf/bin_size_in_bits.hpp>
#include <raptor/build/index_factory.hpp>
#include <raptor/build/max_count_per_partition.hpp>
#include <raptor/build/partition_config.hpp>
#include <raptor/build/store_index.hpp>

namespace raptor
{

void build_ibf(build_arguments const & arguments)
{
    if (arguments.parts == 1u)
    {
        index_factory factory{arguments};
        auto index = factory();
        store_index(arguments.out_path, std::move(index), arguments);
    }
    else
    {
        partition_config const cfg{arguments.parts};
        std::vector<size_t> const kmers_per_partition = max_count_per_partition(cfg, arguments);
        build_arguments args = arguments;
        index_factory factory{args, cfg};

        for (size_t part = 0; part < arguments.parts; ++part)
        {
            args.bits = hibf::bin_size_in_bits(args, kmers_per_partition[part]);
            auto index = factory(part);
            std::filesystem::path out_path{arguments.out_path};
            out_path += "_" + std::to_string(part);
            store_index(out_path, std::move(index), arguments);
        }
    }
}

} // namespace raptor

// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <seqan3/search/views/minimiser_hash.hpp>

#include <raptor/adjust_seed.hpp>
#include <raptor/build/hibf/compute_kmers.hpp>
#include <raptor/dna4_traits.hpp>

namespace raptor::hibf
{

void compute_kmers(robin_hood::unordered_flat_set<size_t> & kmers,
                   build_arguments const & arguments,
                   chopper_pack_record const & record)
{
    if (arguments.is_minimiser)
    {
        uint64_t minimiser_value{};
        for (auto const & filename : record.filenames)
        {
            std::ifstream infile{filename, std::ios::binary};

            while (infile.read(reinterpret_cast<char *>(&minimiser_value), sizeof(minimiser_value)))
                kmers.insert(minimiser_value);
        }
    }
    else
    {
        using sequence_file_t = seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>>;
        for (auto const & filename : record.filenames)
            for (auto && [seq] : sequence_file_t{filename})
                for (auto hash :
                     seq
                         | seqan3::views::minimiser_hash(arguments.shape,
                                                         seqan3::window_size{arguments.window_size},
                                                         seqan3::seed{adjust_seed(arguments.shape.count())}))
                    kmers.insert(hash);
    }
}

} // namespace raptor::hibf

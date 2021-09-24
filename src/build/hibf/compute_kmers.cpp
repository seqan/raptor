// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <seqan3/search/views/minimiser_hash.hpp>

#include <raptor/adjust_seed.hpp>
#include <raptor/build/hibf/compute_kmers.hpp>
#include <raptor/dna4_traits.hpp>

namespace raptor::hibf
{

robin_hood::unordered_flat_set<size_t> compute_kmers(build_config const & config, chopper_pack_record const & record)
{
    robin_hood::unordered_flat_set<size_t> kmers{};

    using sequence_file_t = seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>>;

    for (auto const & filename : record.filenames)
        for (auto && [seq] : sequence_file_t{filename})
            for (auto hash : seq | seqan3::views::minimiser_hash(seqan3::ungapped{config.k},
                                                                 seqan3::window_size{config.k},
                                                                 seqan3::seed{adjust_seed(config.k)}))
                kmers.insert(hash);

    return kmers;
}

void compute_kmers(robin_hood::unordered_flat_set<size_t> & kmers,
                   build_config const & config,
                   chopper_pack_record const & record)
{

    using sequence_file_t = seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>>;
    for (auto const & filename : record.filenames)
        for (auto && [seq] : sequence_file_t{filename})
            for (auto hash : seq | seqan3::views::minimiser_hash(seqan3::ungapped{config.k},
                                                                 seqan3::window_size{config.k},
                                                                 seqan3::seed{adjust_seed(config.k)}))
                kmers.insert(hash);
}

} // namespace raptor::hibf

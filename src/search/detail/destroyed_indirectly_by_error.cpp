// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <random>

#include <seqan3/alphabet/nucleotide/dna4.hpp>

#include <raptor/search/detail/destroyed_indirectly_by_error.hpp>
#include <raptor/search/detail/forward_strand_minimiser.hpp>

namespace raptor::detail
{

std::vector<double> destroyed_indirectly_by_error(size_t const pattern_size,
                                                  size_t const window_size,
                                                  seqan3::shape const shape)
{
    uint8_t const kmer_size{shape.size()};

    using alphabet_t = seqan3::dna4;
    using rank_type = decltype(seqan3::to_rank(alphabet_t{}));
    rank_type max_rank = seqan3::alphabet_size<alphabet_t> - 1;

    std::mt19937_64 gen(0x1D2B8284D988C4D0);
    std::uniform_int_distribution<> dis(0, max_rank);
    std::uniform_int_distribution<> dis2(0, pattern_size - 1);
    std::vector<uint8_t> mins(pattern_size, false);
    std::vector<uint8_t> minse(pattern_size, false);
    std::vector<double> result(window_size - kmer_size + 1, 0);
    std::vector<alphabet_t> sequence;
    sequence.reserve(pattern_size);

    for (size_t iteration = 0; iteration < 10'000; ++iteration)
    {
        sequence.clear();
        std::fill(mins.begin(), mins.end(), false);
        std::fill(minse.begin(), minse.end(), false);

        for (size_t i = 0; i < pattern_size; ++i)
            sequence.push_back(seqan3::assign_rank_to(dis(gen), alphabet_t{}));

        forward_strand_minimiser mini{window{static_cast<uint32_t>(window_size)}, shape};
        mini.compute(sequence);
        for (auto x : mini.minimiser_begin)
            mins[x] = true;

        size_t error_pos = dis2(gen) % pattern_size;
        rank_type new_base = dis(gen) % seqan3::alphabet_size<alphabet_t>;
        while (new_base == seqan3::to_rank(sequence[error_pos]))
            new_base =  dis(gen) % seqan3::alphabet_size<alphabet_t>;
        sequence[error_pos] = seqan3::assign_rank_to(new_base, alphabet_t{});

        mini.compute(sequence);
        for (auto x : mini.minimiser_begin)
            minse[x] = true;

        size_t count = 0;

        for (size_t i = 0; i < pattern_size; ++i)
            count += (mins[i] != minse[i]) && (error_pos < i || i + kmer_size < error_pos);

        ++result[count];
    }

    for (auto & x : result)
        x /= 10'000;

    return result;
}

} // namespace raptor::detail

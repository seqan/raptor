// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <algorithm>
#include <random>

#include <seqan3/search/views/kmer_hash.hpp>

#include <raptor/threshold/forward_strand_minimiser.hpp>

// Begin positions of the minimal (k-mer hash ^ seed) in each window, leftmost on ties, consecutive duplicates removed.
static std::vector<uint64_t>
brute_force(std::vector<seqan3::dna4> const & sequence, uint32_t const window_size, seqan3::shape const shape)
{
    uint64_t const seed{raptor::adjust_seed(shape.count())};
    size_t const kmers_per_window{window_size - shape.size() + 1u};

    std::vector<uint64_t> hashes;
    for (uint64_t const hash : sequence | seqan3::views::kmer_hash(shape))
        hashes.push_back(hash ^ seed);

    std::vector<uint64_t> result;
    for (auto window_begin = hashes.begin(); window_begin + kmers_per_window <= hashes.end(); ++window_begin)
    {
        uint64_t const position =
            std::ranges::min_element(window_begin, window_begin + kmers_per_window) - hashes.begin();
        if (result.empty() || result.back() != position)
            result.push_back(position);
    }
    return result;
}

TEST(forward_strand_minimiser, brute_force)
{
    struct parameters
    {
        size_t query_length;
        uint32_t window_size;
        seqan3::shape shape;
    };

    // Includes w = k, |Q| = w, w > 2k, and a gapped shape.
    std::vector<parameters> const configurations{{250u, 23u, seqan3::ungapped{19u}},
                                                 {250u, 32u, seqan3::ungapped{19u}},
                                                 {100u, 40u, seqan3::ungapped{20u}},
                                                 {50u, 15u, seqan3::ungapped{13u}},
                                                 {15u, 6u, seqan3::ungapped{4u}},
                                                 {20u, 20u, seqan3::ungapped{20u}},
                                                 {65u, 64u, seqan3::ungapped{19u}},
                                                 {64u, 64u, seqan3::ungapped{19u}},
                                                 {100u, 50u, seqan3::ungapped{19u}},
                                                 {100u, 24u, seqan3::shape{seqan3::bin_literal{0b11011101111}}}};

    // The output of the engine is fully specified, unlike the distributions.
    std::mt19937_64 generator{0x5EED};

    for (auto const & [query_length, window_size, shape] : configurations)
    {
        raptor::threshold::forward_strand_minimiser minimiser{raptor::window{window_size}, shape};
        std::vector<seqan3::dna4> sequence(query_length);

        for (size_t run = 0; run < 100u; ++run)
        {
            for (seqan3::dna4 & base : sequence)
                base.assign_rank(static_cast<uint8_t>(generator() % 4u));

            minimiser.compute(sequence);
            EXPECT_EQ(minimiser.minimiser_begin, brute_force(sequence, window_size, shape))
                << "|Q| = " << query_length << ", w = " << window_size << ", shape = " << shape.to_ulong()
                << ", run = " << run;
        }
    }
}

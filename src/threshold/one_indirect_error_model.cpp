// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Implements raptor::threshold::one_indirect_error_model.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <random>

#include <raptor/threshold/forward_strand_minimiser.hpp>
#include <raptor/threshold/one_indirect_error_model.hpp>

namespace raptor::threshold
{

[[nodiscard]] std::vector<double>
one_indirect_error_model(size_t const query_length, size_t const window_size, seqan3::shape const shape)
{
    uint8_t const kmer_size{shape.size()};
    size_t const max_number_of_minimiser{query_length - window_size + 1};
    size_t const iterations{10'000};

    std::mt19937_64 gen{0x1D2B8284D988C4D0};
    std::uniform_int_distribution<size_t> random_error_position{0u, query_length - 1u};
    std::uniform_int_distribution<uint8_t> random_dna4_rank{0u, 3u};
    auto random_dna = [&random_dna4_rank, &gen]()
    {
        return seqan3::assign_rank_to(random_dna4_rank(gen), seqan3::dna4{});
    };

    std::vector<seqan3::dna4> sequence(query_length);

    // Minimiser begin positions of original sequence
    std::vector<uint8_t> minimiser_positions(max_number_of_minimiser, false);
    // Minimiser begin positions after introducing one error into the sequence
    std::vector<uint8_t> minimiser_positions_error(max_number_of_minimiser, false);
    // In the worst case, one error can indirectly affect w minimisers
    std::vector<double> result(window_size + 1, 0.0);
    forward_strand_minimiser fwd_minimiser{window{static_cast<uint32_t>(window_size)}, shape};

    for (size_t iteration = 0; iteration < iterations; ++iteration)
    {
        std::ranges::fill(minimiser_positions, 0u);
        std::ranges::fill(minimiser_positions_error, 0u);
        std::ranges::generate(sequence, random_dna);

        // Minimiser begin positions of original sequence
        fwd_minimiser.compute(sequence);
        for (auto pos : fwd_minimiser.minimiser_begin)
            minimiser_positions[pos] = true;

        // Introduce one error
        size_t const error_position = random_error_position(gen);
        uint8_t new_rank{random_dna4_rank(gen)};
        while (new_rank == seqan3::to_rank(sequence[error_position]))
            new_rank = random_dna4_rank(gen);
        sequence[error_position] = seqan3::assign_rank_to(new_rank, seqan3::dna4{});

        // Minimiser begin positions after introducing one error into the sequence
        fwd_minimiser.compute(sequence);
        for (auto pos : fwd_minimiser.minimiser_begin)
            minimiser_positions_error[pos] = true;

        // Determine number of affected minimisers
        size_t affected_minimiser{};
        // An error destroyed a minimiser indirectly iff
        // (1) A minimiser begin position changed and
        // (2) The error occurs before the window or after the window
        for (size_t i = 0; i < max_number_of_minimiser; ++i)
        {
            affected_minimiser += (minimiser_positions[i] != minimiser_positions_error[i]) && // (1)
                                  ((error_position < i) || (i + kmer_size < error_position)); // (2)
        }

        ++result[affected_minimiser];
    }

    // Convert counts to log probabilities
    double const log_iterations{std::log(iterations)};
    for (double & x : result)
        x = std::log(x) - log_iterations;

    return result;
}

} // namespace raptor::threshold

// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides raptor::threshold::forward_strand_minimiser.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <deque>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

#include <raptor/adjust_seed.hpp>
#include <raptor/strong_types.hpp>

namespace raptor::threshold
{

// Minimiser without looking at reverse complement
struct forward_strand_minimiser
{
private:
    //!\brief The window size of the minimiser.
    uint64_t window_size{};
    //!\brief The shape to use.
    seqan3::shape shape{};
    //!\brief The size of the shape.
    uint8_t shape_size{};
    //!\brief Random but fixed value to xor k-mers with. Counteracts consecutive minimisers.
    uint64_t seed{};
    //!\brief Stores the k-mer hashes of the forward strand.
    std::vector<uint64_t> forward_hashes{};

public:
    //!\brief Stores the begin positions of the minimisers.
    std::vector<uint64_t> minimiser_begin;

    forward_strand_minimiser() = default;                                             //!< Defaulted
    forward_strand_minimiser(forward_strand_minimiser const &) = default;             //!< Defaulted
    forward_strand_minimiser(forward_strand_minimiser &&) = default;                  //!< Defaulted
    forward_strand_minimiser & operator=(forward_strand_minimiser const &) = default; //!< Defaulted
    forward_strand_minimiser & operator=(forward_strand_minimiser &&) = default;      //!< Defaulted
    ~forward_strand_minimiser() = default;                                            //!< Defaulted

    /*!\brief Constructs a minimiser from given k-mer, window size and a seed.
     * \param[in] window_size_ The window size.
     * \param[in] shape_       The shape.
     */
    forward_strand_minimiser(window const window_size_, seqan3::shape const shape_) :
        window_size{window_size_.v},
        shape{shape_},
        shape_size{shape.size()},
        seed{adjust_seed(shape.count())}
    {
        assert(window_size >= shape_size);
    }

    /*!\brief Resize the minimiser.
     * \param[in] window_size_ The new window size.
     * \param[in] shape_       The new shape.
     */
    void resize(window const window_size_, seqan3::shape const shape_)
    {
        window_size = window_size_.v;
        shape = shape_;
        shape_size = shape.size();
        seed = adjust_seed(shape.count());
        assert(window_size >= shape_size);
    }

    void compute(std::vector<seqan3::dna4> const & text)
    {
        assert(window_size && shape_size && seed); // Forgot to initialise/resize?

        size_t const text_length = text.size();
        assert(shape_size <= text_length);
        assert(window_size <= text_length);

        uint64_t const max_number_of_minimiser = text_length - window_size + 1u;
        uint64_t const kmers_per_window = window_size - shape_size + 1u;

        minimiser_begin.clear();
        minimiser_begin.reserve(max_number_of_minimiser);

        // Compute all k-mer hashes.
        auto apply_xor = [this](uint64_t const value)
        {
            return value ^ seed;
        };
        auto kmer_view = text | seqan3::views::kmer_hash(shape) | std::views::transform(apply_xor);
        forward_hashes.assign(kmer_view.begin(), kmer_view.end());

        // Stores hash and begin for all k-mers in the window
        std::deque<std::pair<uint64_t, uint64_t>> window_hashes;

        // Initialisation. We need to compute all hashes for the first window.
        for (uint64_t i = 0; i < kmers_per_window; ++i)
            window_hashes.emplace_back(forward_hashes[i], i);

        // The minimum hash is the minimiser. Store the begin position.
        auto min = std::min_element(std::begin(window_hashes), std::end(window_hashes));
        minimiser_begin.push_back(min->second);

        // For the following windows, we remove the first window k-mer (is now not in window) and add the new k-mer
        // that results from the window shifting.
        for (uint64_t i = kmers_per_window; i < max_number_of_minimiser; ++i)
        {
            // Already store the new hash without removing the first one.
            uint64_t const new_hash{forward_hashes[i + kmers_per_window - 1]}; // Already did kmers_per_window - 1 many
            window_hashes.emplace_back(new_hash, i);

            if (new_hash < min->second) // New hash is the minimum.
            {
                min = std::prev(std::end(window_hashes));
                minimiser_begin.push_back(min->second);
            }
            else if (min == std::begin(window_hashes)) // Minimum is the yet-to-be-removed begin of the window.
            {
                // The first hash will be removed, the last one is caught by the previous if.
                min = std::min_element(++std::begin(window_hashes), std::prev(std::end(window_hashes)));
                minimiser_begin.push_back(min->second);
            }

            window_hashes.pop_front(); // Remove the first k-mer.
        }
        return;
    }
};

} // namespace raptor::threshold

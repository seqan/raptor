// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

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

        struct kmer
        {
            uint64_t hash{};
            uint64_t position{};

            constexpr auto operator<=>(kmer const & other) const = default;
        };

        // Stores hash and begin for all k-mers in the window
        std::deque<kmer> window_hashes;

        // Initialisation. We need to compute all hashes for the first window.
        for (uint64_t i = 0; i < kmers_per_window; ++i)
            window_hashes.emplace_back(kmer{.hash = forward_hashes[i], .position = i});

        // The minimum hash is the minimiser. Store the begin position.
        auto min = std::ranges::min(window_hashes);
        minimiser_begin.push_back(min.position);

        // For the following windows, we remove the first window k-mer (is now not in window) and add the new k-mer
        // that results from the window shifting.
        for (uint64_t i = kmers_per_window; i < max_number_of_minimiser; ++i)
        {
            uint64_t const new_hash{forward_hashes[i + kmers_per_window - 1]}; // Already did kmers_per_window - 1 many

            // There are two conditions when we need to recompute the minimum:
            bool const minimiser_leaves_window = window_hashes.front() == min;
            bool const new_hash_is_min = new_hash < min.hash;

            // Rolling hash / sliding window
            window_hashes.pop_front();
            window_hashes.emplace_back(kmer{.hash = new_hash, .position = i});

            // Update the minimum
            if (new_hash_is_min || minimiser_leaves_window)
            {
                min = new_hash_is_min ? window_hashes.back() : std::ranges::min(window_hashes);
                minimiser_begin.push_back(min.position);
            }
        }
        return;
    }
};

} // namespace raptor::threshold

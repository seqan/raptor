#pragma once

#include <numeric>
#include <random>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

#include <raptor/shared.hpp>

namespace raptor::detail
{

// Minimiser without looking at reverse complement
struct forward_strand_minimiser
{
private:
    //!\brief The alphabet type.
    using alphabet_t = seqan3::dna4;
    //!\brief The text type.
    using text_t = seqan3::dna4_vector;

    //!\brief The window size of the minimiser.
    uint64_t w{26};
    //!\brief The size of the k-mers.
    uint8_t k{20};
    //!\brief Random but fixed value to xor k-mers with. Counteracts consecutive minimisers.
    uint64_t seed{adjust_seed(k)};

    //!\brief Stores the k-mer hashes of the forward strand.
    std::vector<uint64_t> forward_hashes;

public:

    //!\brief Stores the hashes of the minimisers.
    std::vector<uint64_t> minimiser_hash;
    //!\brief Stores the begin positions of the minimisers.
    std::vector<uint64_t> minimiser_begin;
    //!\brief Stores the end positions of the minimisers.
    std::vector<uint64_t> minimiser_end;

    forward_strand_minimiser() = default;                                        //!< Defaulted
    forward_strand_minimiser(forward_strand_minimiser const &) = default;                       //!< Defaulted
    forward_strand_minimiser(forward_strand_minimiser &&) = default;                            //!< Defaulted
    forward_strand_minimiser & operator=(forward_strand_minimiser const &) = default;           //!< Defaulted
    forward_strand_minimiser & operator=(forward_strand_minimiser &&) = default;                //!< Defaulted
    ~forward_strand_minimiser() = default;                                       //!< Defaulted

    /*!\brief Constructs a minimiser from given k-mer, window size and a seed.
     * \param[in] w_    The window size.
     * \param[in] k_    The k-mer size.
     * \param[in] seed_ The seed to use. Default: 0x8F3F73B5CF1C9ADE.
     */
    forward_strand_minimiser(window const w_, kmer const k_, uint64_t const seed_ = 0x8F3F73B5CF1C9ADE) :
        w{w_.v}, k{k_.v}, seed{adjust_seed(k_.v, seed_)}
    {}

    /*!\brief Resize the minimiser.
     * \param[in] w_    The new window size.
     * \param[in] k_    The new k-mer size.
     * \param[in] seed_ The new seed to use. Default: 0x8F3F73B5CF1C9ADE.
     */
    void resize(window const w_, kmer const k_, uint64_t const seed_ = 0x8F3F73B5CF1C9ADE)
    {
        w = w_.v;
        k = k_.v;
        seed = adjust_seed(k_.v, seed_);
    }

    void compute(text_t const & text)
    {
        uint64_t text_length = std::ranges::size(text);

        forward_hashes.clear();
        minimiser_hash.clear();
        minimiser_begin.clear();
        minimiser_end.clear();

        // Return empty vector if text is shorter than k.
        if (k > text_length)
            return;

        uint64_t possible_minimisers = text_length > w ? text_length - w + 1u : 1u;
        assert(w >= k);
        uint64_t kmers_per_window = w - k + 1u;

        // Helper lambda for xor'ing values depending on `do_xor`.
        auto apply_xor = [this] (uint64_t const val)
        {
            return val ^ seed;
        };

        // Compute all k-mer hashes for both forward and reverse strand.
        forward_hashes = text |
                         seqan3::views::kmer_hash(seqan3::ungapped{k}) |
                         std::views::transform(apply_xor) |
                         seqan3::views::to<std::vector<uint64_t>>;

        // Choose the minimisers.
        minimiser_hash.reserve(possible_minimisers);
        minimiser_begin.reserve(possible_minimisers);
        minimiser_end.reserve(possible_minimisers);

        // Stores hash, begin and end for all k-mers in the window
        std::deque<std::tuple<uint64_t, uint64_t, uint64_t>> window_values;

        // Initialisation. We need to compute all hashes for the first window.
        for (uint64_t i = 0; i < kmers_per_window; ++i)
            window_values.emplace_back(forward_hashes[i], i, i + k - 1);

        auto min = std::min_element(std::begin(window_values), std::end(window_values));
        minimiser_hash.push_back(std::get<0>(*min));
        minimiser_begin.push_back(std::get<1>(*min));
        minimiser_end.push_back(std::get<2>(*min));

        // For the following windows, we remove the first window k-mer (is now not in window) and add the new k-mer
        // that results from the window shifting
        bool minimiser_changed{false};
        for (uint64_t i = 1; i < possible_minimisers; ++i)
        {
            // Shift the window.
            // If current minimiser leaves the window, we need to decide on a new one.
            if (min == std::begin(window_values))
            {
                window_values.pop_front();
                min = std::min_element(std::begin(window_values), std::end(window_values));
                minimiser_changed = true;
            }
            else
            {
                window_values.pop_front();
            }

            window_values.emplace_back(forward_hashes[kmers_per_window - 1 + i],
                                       kmers_per_window + i - 1,
                                       kmers_per_window + i + k - 2);

            if (std::get<0>(window_values.back()) < std::get<0>(*min))
            {
                min = std::prev(std::end(window_values));
                minimiser_changed = true;
            }

            if (minimiser_changed)
            {
                minimiser_hash.push_back(std::get<0>(*min));
                minimiser_begin.push_back(std::get<1>(*min));
                minimiser_end.push_back(std::get<2>(*min));
                minimiser_changed = false;
            }
        }
        return;
    }
};

} // namespace raptor::detail


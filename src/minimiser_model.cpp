// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/eseiler/minimizer_thresholds/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides stuff.
 */

#include <random>
#include <vector>

#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/range/views/complement.hpp>
#include <seqan3/range/views/kmer_hash.hpp>

#if SEQAN3_WITH_CEREAL
#include <cereal/types/vector.hpp>
#else
#error "Cereal not found, make sure you cloned the project recursively."
#endif

#include <shared.hpp>

struct minimizer
{
private:
    //!\brief The alphabet type.
    using alphabet_t = seqan3::dna4;
    //!\brief The text type.
    using text_t = seqan3::dna4_vector;

    //!\brief The window size of the minimizer.
    uint64_t w{26};
    //!\brief The size of the k-mers.
    uint8_t k{20};
    //!\brief Random but fixed value to xor k-mers with. Counteracts consecutive minimizers.
    uint64_t seed{0x8F3F73B5CF1C9ADE};

    //!\brief Stores the k-mer hashes of the forward strand.
    std::vector<uint64_t> forward_hashes;
    //!\brief Stores the k-mer hashes of the reverse complement strand.
    std::vector<uint64_t> reverse_hashes;


public:

    //!\brief Stores the hashes of the minimizers.
    std::vector<uint64_t> minimizer_hash;
    //!\brief Stores the begin positions of the minimizers.
    std::vector<uint64_t> minimizer_begin;
    //!\brief Stores the end positions of the minimizers.
    std::vector<uint64_t> minimizer_end;

    minimizer() = default;                                        //!< Defaulted
    minimizer(minimizer const &) = default;                       //!< Defaulted
    minimizer(minimizer &&) = default;                            //!< Defaulted
    minimizer & operator=(minimizer const &) = default;           //!< Defaulted
    minimizer & operator=(minimizer &&) = default;                //!< Defaulted
    ~minimizer() = default;                                       //!< Defaulted

    /*!\brief Constructs a minimizer from given k-mer, window size and a seed.
     * \param[in] w_    The window size.
     * \param[in] k_    The k-mer size.
     * \param[in] seed_ The seed to use. Default: 0x8F3F73B5CF1C9ADE.
     */
    minimizer(window const w_, kmer const k_, uint64_t const seed_ = 0x8F3F73B5CF1C9ADE) :
        w{w_.v}, k{k_.v}, seed{seed_}
    {}

    /*!\brief Resize the minimizer.
     * \param[in] w_    The new window size.
     * \param[in] k_    The new k-mer size.
     * \param[in] seed_ The new seed to use. Default: 0x8F3F73B5CF1C9ADE.
     */
    void resize(window const w_, kmer const k_, uint64_t const seed_ = 0x8F3F73B5CF1C9ADE)
    {
        w = w_.v;
        k = k_.v;
        seed = seed_;
    }

    std::vector<uint64_t> hashes(text_t const & text)
    {
        compute(text);
        return minimizer_hash;
    }

    std::vector<uint64_t> hashes_multi(text_t const & text)
    {
        compute_multi(text);
        return minimizer_hash;
    }

    void compute(text_t const & text)
    {
        uint64_t text_length = std::ranges::size(text);

        forward_hashes.clear();
        reverse_hashes.clear();
        minimizer_hash.clear();
        minimizer_begin.clear();
        minimizer_end.clear();

        // Return empty vector if text is shorter than k.
        if (k > text_length)
            return;

        auto rc_text = text | seqan3::views::complement | std::views::reverse;
        uint64_t possible_minimizers = text_length > w ? text_length - w + 1u : 1u;
        uint64_t possible_kmers = text_length - k + 1;
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

        reverse_hashes = rc_text |
                         seqan3::views::kmer_hash(seqan3::ungapped{k}) |
                         std::views::transform(apply_xor) |
                         seqan3::views::to<std::vector<uint64_t>>;

        // Choose the minimizers.
        minimizer_hash.reserve(possible_minimizers);
        minimizer_begin.reserve(possible_minimizers);
        minimizer_end.reserve(possible_minimizers);

        // Stores hash, begin and end for all k-mers in the window
        std::deque<std::tuple<uint64_t, uint64_t, uint64_t>> window_values;

        // Initialisation. We need to compute all hashes for the first window.
        for (uint64_t i = 0; i < kmers_per_window; ++i)
        {
            // Get smallest canonical k-mer.
            uint64_t forward_hash = forward_hashes[i];
            uint64_t reverse_hash = reverse_hashes[possible_kmers - i - 1];
            window_values.emplace_back(std::min(forward_hash, reverse_hash), i, i + k - 1);
        }

        auto min = std::min_element(std::begin(window_values), std::end(window_values));
        minimizer_hash.push_back(std::get<0>(*min));
        minimizer_begin.push_back(std::get<1>(*min));
        minimizer_end.push_back(std::get<2>(*min));


        // For the following windows, we remove the first window k-mer (is now not in window) and add the new k-mer
        // that results from the window shifting
        bool minimizer_changed{false};
        for (uint64_t i = 1; i < possible_minimizers; ++i)
        {
            // Shift the window.
            // If current minimizer leaves the window, we need to decide on a new one.
            if (min == std::begin(window_values))
            {
                window_values.pop_front();
                min = std::min_element(std::begin(window_values), std::end(window_values));
                minimizer_changed = true;
            }
            else
            {
                window_values.pop_front();
            }

            uint64_t forward_hash = forward_hashes[kmers_per_window - 1 + i];
            uint64_t reverse_hash = reverse_hashes[possible_kmers - kmers_per_window - i];
            window_values.emplace_back(std::min(forward_hash, reverse_hash),
                                       kmers_per_window + i - 1,
                                       kmers_per_window + i + k - 2);

            if (std::get<0>(window_values.back()) < std::get<0>(*min))
            {
                min = std::prev(std::end(window_values));
                minimizer_changed = true;
            }

            if (minimizer_changed)
            {
                minimizer_hash.push_back(std::get<0>(*min));
                minimizer_begin.push_back(std::get<1>(*min));
                minimizer_end.push_back(std::get<2>(*min));
                minimizer_changed = false;
            }
        }
        return;
    }

    void compute_multi(text_t const & text)
    {
        uint64_t text_length = std::ranges::size(text);

        forward_hashes.clear();
        reverse_hashes.clear();
        minimizer_hash.clear();
        minimizer_begin.clear();
        minimizer_end.clear();

        // Return empty vector if text is shorter than k.
        if (k > text_length)
            return;

        auto rc_text = text | seqan3::views::complement | std::views::reverse;
        uint64_t possible_minimizers = text_length > w ? text_length - w + 1u : 1u;
        uint64_t possible_kmers = text_length - k + 1;
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

        reverse_hashes = rc_text |
                         seqan3::views::kmer_hash(seqan3::ungapped{k}) |
                         std::views::transform(apply_xor) |
                         seqan3::views::to<std::vector<uint64_t>>;

        // Choose the minimizers.
        minimizer_hash.reserve(possible_minimizers);

        // Stores hash, begin and end for all k-mers in the window
        std::deque<std::tuple<uint64_t, uint64_t, uint64_t>> window_values;

        // Initialisation. We need to compute all hashes for the first window.
        for (uint64_t i = 0; i < kmers_per_window; ++i)
        {
            // Get smallest canonical k-mer.
            uint64_t forward_hash = forward_hashes[i];
            uint64_t reverse_hash = reverse_hashes[possible_kmers - i - 1];
            window_values.emplace_back(std::min(forward_hash, reverse_hash), i, i + k - 1);
        }

        auto min = std::min_element(std::begin(window_values), std::end(window_values));
        minimizer_hash.push_back(std::get<0>(*min));

        // For the following windows, we remove the first window k-mer (is now not in window) and add the new k-mer
        // that results from the window shifting
        for (uint64_t i = 1; i < possible_minimizers; ++i)
        {
            // Shift the window.
            // If current minimizer leaves the window, we need to decide on a new one.
            if (min == std::begin(window_values))
            {
                window_values.pop_front();
                min = std::min_element(std::begin(window_values), std::end(window_values));
            }
            else
            {
                window_values.pop_front();
            }

            uint64_t forward_hash = forward_hashes[kmers_per_window - 1 + i];
            uint64_t reverse_hash = reverse_hashes[possible_kmers - kmers_per_window - i];
            window_values.emplace_back(std::min(forward_hash, reverse_hash),
                                       kmers_per_window + i - 1,
                                       kmers_per_window + i + k - 2);

            if (std::get<0>(window_values.back()) < std::get<0>(*min))
            {
                min = std::prev(std::end(window_values));
            }

            minimizer_hash.push_back(std::get<0>(*min));
        }
        return;
    }

};

std::vector<size_t> pascal_row(size_t n)
{
    std::vector<size_t> result(n + 1);
    result[0] = 1;

    for (size_t i = 1; i <= n; ++i)
        result[i] = result[i - 1] * (n + 1 - i) / i;

    return result;
}

std::tuple<double, std::vector<double>>
simple_model(size_t const kmer_size,
             std::vector<double> const & proba_x,
             std::vector<double> const & indirect_errors)
{
    // Find worst case.
    double max = 0;

    for (size_t i = 0; i < std::ranges::size(proba_x); ++i)
    {
        // Sum up kmer_size many positions.
        double tmp = std::accumulate(proba_x.begin() + i,
                                     proba_x.begin() + std::min(std::ranges::size(proba_x), i + kmer_size),
                                     0.0);

        max = std::max(tmp, max);
    }

    std::vector<size_t> coefficients{pascal_row(kmer_size)};
    std::vector<double> probabilities(kmer_size + 1);
    double p_mean = max / static_cast<double>(kmer_size);
    double p_sum = 0;

    for (size_t i = 0; i <= kmer_size; ++i)
    {
        double p_i_error = coefficients[i] * std::pow(p_mean, i) * std::pow(1 - p_mean, kmer_size - i);

        for (size_t j = 0; j < indirect_errors.size() && i + j <= kmer_size; ++j)
            probabilities[i + j] += p_i_error * indirect_errors[j];

        p_sum += probabilities[i];
    }

    for (auto & x : probabilities)
        x /= p_sum;

    return {p_mean, probabilities};
}

void impl(size_t const minimizers_left,
          std::vector<double> const & proba,
          std::vector<size_t> error_distribution,
          size_t const current_error_index,
          double & result)
{
    if (!minimizers_left)
    {
        double tmp = 1;

        // Probabilities that error_distribution[i] many minimizers are affected by error i.
        for (size_t i = 0; i < current_error_index; ++i)
            tmp *= proba[error_distribution[i]];
        // Then the other errors must not affect any minimizers.
        for (size_t i = current_error_index; i < error_distribution.size(); ++i)
            tmp *= proba[0];

        result += tmp;
        return;
    }

    if (current_error_index >= error_distribution.size())
        return;

    // Enumerate. Can't use too many minimizers and can't destroy more than proba.size() with one error.
    for (size_t i = 0; i <= minimizers_left && i < proba.size(); ++i)
    {
        error_distribution[current_error_index] = i;
        impl(minimizers_left - i, proba, error_distribution, current_error_index + 1, result);
    }
}

double enumerate_all_errors(size_t const number_of_minimizers, size_t const errors, std::vector<double> const & proba)
{
    double result = 0;
    impl(number_of_minimizers, proba, std::vector<size_t>(errors, 0), 0, result);
    return result;
}

// Minimizer without looking at reverse complement
struct boring_minimizer
{
private:
    //!\brief The alphabet type.
    using alphabet_t = seqan3::dna4;
    //!\brief The text type.
    using text_t = seqan3::dna4_vector;

    //!\brief The window size of the minimizer.
    uint64_t w{26};
    //!\brief The size of the k-mers.
    uint8_t k{20};
    //!\brief Random but fixed value to xor k-mers with. Counteracts consecutive minimizers.
    uint64_t seed{adjust_seed(k)};

    //!\brief Stores the k-mer hashes of the forward strand.
    std::vector<uint64_t> forward_hashes;

public:

    //!\brief Stores the hashes of the minimizers.
    std::vector<uint64_t> minimizer_hash;
    //!\brief Stores the begin positions of the minimizers.
    std::vector<uint64_t> minimizer_begin;
    //!\brief Stores the end positions of the minimizers.
    std::vector<uint64_t> minimizer_end;

    boring_minimizer() = default;                                        //!< Defaulted
    boring_minimizer(boring_minimizer const &) = default;                       //!< Defaulted
    boring_minimizer(boring_minimizer &&) = default;                            //!< Defaulted
    boring_minimizer & operator=(boring_minimizer const &) = default;           //!< Defaulted
    boring_minimizer & operator=(boring_minimizer &&) = default;                //!< Defaulted
    ~boring_minimizer() = default;                                       //!< Defaulted

    /*!\brief Constructs a minimizer from given k-mer, window size and a seed.
     * \param[in] w_    The window size.
     * \param[in] k_    The k-mer size.
     * \param[in] seed_ The seed to use. Default: 0x8F3F73B5CF1C9ADE.
     */
    boring_minimizer(window const w_, kmer const k_, uint64_t const seed_ = 0x8F3F73B5CF1C9ADE) :
        w{w_.v}, k{k_.v}, seed{adjust_seed(k_.v, seed_)}
    {}

    /*!\brief Resize the minimizer.
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
        minimizer_hash.clear();
        minimizer_begin.clear();
        minimizer_end.clear();

        // Return empty vector if text is shorter than k.
        if (k > text_length)
            return;

        uint64_t possible_minimizers = text_length > w ? text_length - w + 1u : 1u;
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

        // Choose the minimizers.
        minimizer_hash.reserve(possible_minimizers);
        minimizer_begin.reserve(possible_minimizers);
        minimizer_end.reserve(possible_minimizers);

        // Stores hash, begin and end for all k-mers in the window
        std::deque<std::tuple<uint64_t, uint64_t, uint64_t>> window_values;

        // Initialisation. We need to compute all hashes for the first window.
        for (uint64_t i = 0; i < kmers_per_window; ++i)
            window_values.emplace_back(forward_hashes[i], i, i + k - 1);

        auto min = std::min_element(std::begin(window_values), std::end(window_values));
        minimizer_hash.push_back(std::get<0>(*min));
        minimizer_begin.push_back(std::get<1>(*min));
        minimizer_end.push_back(std::get<2>(*min));

        // For the following windows, we remove the first window k-mer (is now not in window) and add the new k-mer
        // that results from the window shifting
        bool minimizer_changed{false};
        for (uint64_t i = 1; i < possible_minimizers; ++i)
        {
            // Shift the window.
            // If current minimizer leaves the window, we need to decide on a new one.
            if (min == std::begin(window_values))
            {
                window_values.pop_front();
                min = std::min_element(std::begin(window_values), std::end(window_values));
                minimizer_changed = true;
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
                minimizer_changed = true;
            }

            if (minimizer_changed)
            {
                minimizer_hash.push_back(std::get<0>(*min));
                minimizer_begin.push_back(std::get<1>(*min));
                minimizer_end.push_back(std::get<2>(*min));
                minimizer_changed = false;
            }
        }
        return;
    }
};

template <seqan3::alphabet alphabet_t>
std::vector<double> destroyed_indirectly_by_error(size_t const pattern_size, size_t const window_size, uint8_t const kmer_size)
{
    using rank_type = decltype(seqan3::to_rank(alphabet_t{}));
    rank_type max_rank = seqan3::alphabet_size<alphabet_t> - 1;

    std::mt19937_64 gen(0x1D2B8284D988C4D0);
    std::uniform_int_distribution<> dis(0, max_rank);
    std::uniform_int_distribution<> dis2(0, pattern_size - 1);
    std::vector<uint8_t> mins(pattern_size, false);
    std::vector<uint8_t> minse(pattern_size, false);
    std::vector<double> result(window_size - kmer_size, 0);
    std::vector<alphabet_t> sequence;
    sequence.reserve(pattern_size);

    for (size_t iteration = 0; iteration < 10'000; ++iteration)
    {
        sequence.clear();
        std::fill(mins.begin(), mins.end(), false);
        std::fill(minse.begin(), minse.end(), false);

        for (size_t i = 0; i < pattern_size; ++i)
            sequence.push_back(seqan3::assign_rank_to(dis(gen), alphabet_t{}));

        boring_minimizer mini{window{window_size}, kmer{kmer_size}};
        mini.compute(sequence);
        for (auto x : mini.minimizer_begin)
            mins[x] = true;

        size_t error_pos = dis2(gen) % pattern_size;
        rank_type new_base = dis(gen) % seqan3::alphabet_size<alphabet_t>;
        while (new_base == seqan3::to_rank(sequence[error_pos]))
            new_base =  dis(gen) % seqan3::alphabet_size<alphabet_t>;
        sequence[error_pos] = seqan3::assign_rank_to(new_base, alphabet_t{});

        mini.compute(sequence);
        for (auto x : mini.minimizer_begin)
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

std::vector<size_t> precompute_threshold(size_t const pattern_size,
                                         size_t const window_size,
                                         uint8_t const kmer_size,
                                         size_t const errors,
                                         double const tau)
{
    if (window_size == kmer_size)
        return {pattern_size + 1 > (errors + 1) * kmer_size ? pattern_size + 1 - (errors + 1) * kmer_size : 0};

    std::vector<size_t> thresholds;
    size_t const kmers_per_window = window_size - kmer_size + 1;
    size_t const kmers_per_pattern = pattern_size - kmer_size + 1;

    size_t const minimal_number_of_minimizers = std::ceil(kmers_per_pattern / static_cast<double>(kmers_per_window));
    size_t const maximal_number_of_minimizers = pattern_size - window_size + 1;

    std::vector<double> indirect_errors;
    indirect_errors = destroyed_indirectly_by_error<seqan3::dna4>(pattern_size, window_size, kmer_size);

    // Iterate over the possible number of minimizers
    for (size_t number_of_minimizers = minimal_number_of_minimizers; number_of_minimizers <= maximal_number_of_minimizers; ++number_of_minimizers)
    {
        std::vector<double> proba_x(kmers_per_pattern, number_of_minimizers / static_cast<double>(kmers_per_pattern));

        auto [p_mean, proba] = simple_model(kmer_size, proba_x, indirect_errors);
        (void) p_mean;

        std::vector<double> proba_error(number_of_minimizers, 0);
        for (size_t i = 0; i < number_of_minimizers; ++i)
            proba_error[i] = enumerate_all_errors(i, errors, proba);

        double sum = std::accumulate(proba_error.begin(), proba_error.end(), 0.0);
        for (auto & x : proba_error)
            x /= sum;

        double n =0;
        for (size_t i = 0; i < number_of_minimizers; ++i)
        {
            n += proba_error[i];

            if (n >= tau)
            {
                thresholds.push_back(number_of_minimizers - i);
                break;
            }
        }
    }
    assert(thresholds.size() != 0);
    return thresholds;
}

void do_cerealisation_out(std::vector<size_t> const & vec, search_arguments const & arguments)
{
    std::filesystem::path filename = arguments.ibf_file.parent_path() / ("binary_p" + std::to_string(arguments.pattern_size) +
                                                                         "_w" + std::to_string(arguments.window_size) +
                                                                         "_k" + std::to_string(arguments.kmer_size) +
                                                                         "_e" + std::to_string(arguments.errors) +
                                                                         "_tau" + std::to_string(arguments.tau));
    std::ofstream os{filename, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(vec);
}

bool do_cerealisation_in(std::vector<size_t> & vec, search_arguments const & arguments)
{
    std::filesystem::path filename = arguments.ibf_file.parent_path() / ("binary_p" + std::to_string(arguments.pattern_size) +
                                                                         "_w" + std::to_string(arguments.window_size) +
                                                                         "_k" + std::to_string(arguments.kmer_size) +
                                                                         "_e" + std::to_string(arguments.errors) +
                                                                         "_tau" + std::to_string(arguments.tau));
    if (!std::filesystem::exists(filename))
        return false;

    std::ifstream is{filename, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    iarchive(vec);
    return true;
}

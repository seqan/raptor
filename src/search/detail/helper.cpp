// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <numeric>
#include <random>

#include <raptor/search/detail/forward_strand_minimiser.hpp>
#include <raptor/search/detail/helper.hpp>

namespace raptor::detail
{

std::vector<size_t> pascal_row(size_t n)
{
    std::vector<size_t> result(n + 1);
    result[0] = 1;

    for (size_t i = 1; i <= n; ++i)
        result[i] = result[i - 1] * (n + 1 - i) / i;

    return result;
}

std::tuple<double, std::vector<double>>
simple_model(size_t const kmer_size, std::vector<double> const & proba_x, std::vector<double> const & indirect_errors)
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

void impl(size_t const minimisers_left,
          std::vector<double> const & proba,
          std::vector<size_t> error_distribution,
          size_t const current_error_index,
          double & result)
{
    if (!minimisers_left)
    {
        double tmp = 1;

        // Probabilities that error_distribution[i] many minimisers are affected by error i.
        for (size_t i = 0; i < current_error_index; ++i)
            tmp *= proba[error_distribution[i]];
        // Then the other errors must not affect any minimisers.
        for (size_t i = current_error_index; i < error_distribution.size(); ++i)
            tmp *= proba[0];

        result += tmp;
        return;
    }

    if (current_error_index >= error_distribution.size())
        return;

    // Enumerate. Can't use too many minimisers and can't destroy more than proba.size() with one error.
    for (size_t i = 0; i <= minimisers_left && i < proba.size(); ++i)
    {
        error_distribution[current_error_index] = i;
        impl(minimisers_left - i, proba, error_distribution, current_error_index + 1, result);
    }
}

double enumerate_all_errors(size_t const number_of_minimisers, size_t const errors, std::vector<double> const & proba)
{
    double result = 0;
    impl(number_of_minimisers, proba, std::vector<size_t>(errors, 0), 0, result);
    return result;
}

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

        forward_strand_minimiser mini{window{window_size}, shape};
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

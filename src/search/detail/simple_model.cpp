// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <cmath>
#include <numeric>
#include <seqan3/std/ranges>

#include <raptor/search/detail/pascal_row.hpp>
#include <raptor/search/detail/simple_model.hpp>

namespace raptor::detail
{

std::tuple<double, std::vector<double>> simple_model(size_t const kmer_size,
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


} // namespace raptor::detail

// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <raptor/threshold/logspace.hpp>
#include <raptor/threshold/pascal_row.hpp>
#include <raptor/threshold/one_error_model.hpp>

namespace raptor::threshold
{

[[nodiscard]] std::vector<double> one_error_model(size_t const kmer_size,
                                                  double const p_mean,
                                                  std::vector<double> const & affected_by_one_error_indirectly_prob)
{
    std::vector<double> const coefficients{pascal_row(kmer_size)};
    // Probabilities that i minimisers are affected by one error.
    std::vector<double> probabilities(kmer_size + 1, logspace::negative_inf);
    double const inv_p_mean{logspace::substract(0, p_mean)};
    double p_sum{logspace::negative_inf};

    // One error can affect between 0 and k minimisers.
    for (size_t i = 0; i <= kmer_size; ++i)
    {
        // Probability that i minimisers are directly modified by one error.
        double const p_i{coefficients[i] + i * p_mean + (kmer_size - i) * inv_p_mean};

        // Incorporate indirect errors:
        // p_i = probability that one error affects i minimisers directly
        // affected_by_one_error_indirectly_prob[j] = probability that one error affects j minimisers indirectly
        // At most k many minimisers can be affected: i + j <= kmer_size
        // indirect and direct errors occur independently
        // The for loops will enumerate all combinations of i and j for achieving 0 to k many affected minimiser.
        for (size_t j = 0; j < affected_by_one_error_indirectly_prob.size() && i + j <= kmer_size; ++j)
            probabilities[i + j] = logspace::add(probabilities[i + j], p_i + affected_by_one_error_indirectly_prob[j]);

        // Keep track of sum of probabilities.
        p_sum = logspace::add(p_sum, probabilities[i]);
    }

    // Normalise probabilities.
    for (double & x : probabilities)
        x -= p_sum;

    return probabilities;
}


} // namespace raptor::threshold

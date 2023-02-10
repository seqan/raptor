// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Implements raptor::threshold::one_error_model.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <cstddef>
#include <numeric>

#include <raptor/threshold/logspace.hpp>
#include <raptor/threshold/one_error_model.hpp>
#include <raptor/threshold/pascal_row.hpp>

namespace raptor::threshold
{

[[nodiscard]] std::vector<double> one_error_model(size_t const kmer_size,
                                                  double const p_mean,
                                                  std::vector<double> const & affected_by_one_error_indirectly_prob)
{
    size_t const window_size{affected_by_one_error_indirectly_prob.size() - 1};
    std::vector<double> const coefficients{pascal_row(kmer_size)};
    // Probabilities that i minimisers are affected by one error.
    std::vector<double> probabilities(window_size + 1, logspace::negative_inf);
    double const inv_p_mean{logspace::substract(0, p_mean)};

    // One error can affect between 0 and w minimisers. Direct errors 0 to k.
    for (size_t i = 0; i <= kmer_size; ++i)
    {
        // Probability that i minimisers are directly modified by one error.
        double const p_direct{coefficients[i] + i * p_mean + (kmer_size - i) * inv_p_mean};

        // Incorporate indirect errors:
        // p_direct = probability that one error affects i minimisers directly
        // affected_by_one_error_indirectly_prob[j] = probability that one error affects j minimisers indirectly
        // At most w many minimisers can be affected: i + j <= window_size
        // indirect and direct errors occur independently
        // The for loops will enumerate all combinations of i and j for achieving 0 to w many affected minimiser.
        for (size_t j = 0; i + j <= window_size; ++j)
            probabilities[i + j] =
                logspace::add(probabilities[i + j], p_direct + affected_by_one_error_indirectly_prob[j]);
    }

    // Normalise probabilities.
    double const p_sum =
        std::accumulate(probabilities.begin(), probabilities.end(), logspace::negative_inf, logspace::add_fn{});

    for (double & x : probabilities)
        x -= p_sum;

    return probabilities;
}

} // namespace raptor::threshold

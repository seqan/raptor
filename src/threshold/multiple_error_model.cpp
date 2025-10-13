// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::threshold::multiple_error_model.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <algorithm> // for clamp
#include <cstddef>   // for size_t
#include <vector>    // for vector

#include <raptor/threshold/logspace.hpp>             // for add, negative_inf
#include <raptor/threshold/multiple_error_model.hpp> // for multiple_error_model

namespace raptor::threshold
{

// affected_by_error [2, 0, 1] => 3 errors. First affects two minimisers, second none, third one minimiser.
// current_error: All possible error configs are enumerated. current_error describes which error is to be enumerated.
void impl(size_t const minimisers_to_affect,
          std::vector<double> const & affected_by_one_error_prob,
          std::vector<size_t> affected_by_error,
          size_t const current_error,
          double & result)
{
    if (!minimisers_to_affect)
    {
        // Determine the probability that this comination of errors occurs.
        double current_prob{};

        // The probability that the errors affect minimisers in this specific way.
        for (size_t i = 0; i < current_error; ++i)
            current_prob += affected_by_one_error_prob[affected_by_error[i]];
        // Then the other errors must not affect any minimisers.
        for (size_t i = current_error; i < affected_by_error.size(); ++i)
            current_prob += affected_by_one_error_prob[0];

        result = logspace::add(result, current_prob);
        return;
    }

    if (current_error >= affected_by_error.size())
        return;

    // Enumerate. Can't use too many minimisers and can't destroy more than proba.size() with one error.
    for (size_t i = 0; i <= minimisers_to_affect && i < affected_by_one_error_prob.size(); ++i)
    {
        affected_by_error[current_error] = i;
        impl(minimisers_to_affect - i, affected_by_one_error_prob, affected_by_error, current_error + 1, result);
    }
}

[[nodiscard]] std::vector<double> multiple_error_model(size_t const number_of_minimisers,
                                                       size_t const errors,
                                                       std::vector<double> const & affected_by_one_error_prob)
{
    size_t const window_size{affected_by_one_error_prob.size() - 1};
    size_t const max_affected{std::clamp<size_t>(errors * window_size, 0u, number_of_minimisers)};
    std::vector<double> affected_by_e_errors(max_affected + 1, 0);

    double sum{logspace::negative_inf};

    // Enumerate all combinations which lead to i many affected minimisers using e errors.
    for (size_t i = 0; i <= max_affected; ++i)
    {
        double result{logspace::negative_inf};
        impl(i, affected_by_one_error_prob, std::vector<size_t>(errors, 0), 0, result);
        affected_by_e_errors[i] = result;
        sum = logspace::add(sum, result);
    }

    // Normalise.
    for (auto & x : affected_by_e_errors)
        x -= sum;

    return affected_by_e_errors;
}

} // namespace raptor::threshold

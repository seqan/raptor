// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <raptor/search/detail/enumerate_all_errors.hpp>

namespace raptor::detail
{

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

} // namespace raptor::detail

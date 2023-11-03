// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides functionality for working with probabilities in logspace.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <cmath>

namespace raptor::logspace
{

constexpr double ln_2{0.693147180559945309417232121458176568L};
constexpr double negative_inf{-std::numeric_limits<double>::infinity()};

//!\brief The log of a sum of two log terms.
// Produces correct result if either term is -inf.
// Needs a check for the case where both terms are -inf.
[[nodiscard]] inline double add(double const log_x, double const log_y) noexcept
{
    double const max{std::max(log_x, log_y)};
    return max == negative_inf ? negative_inf : max + std::log1p(std::exp(-std::abs(log_x - log_y)));
}

//!\brief The log of a sum of multiple log terms.
template <typename... types>
[[nodiscard]] double add(double const log_x, double const log_y, types... logs) noexcept
{
    return add(add(log_y, log_x), logs...);
}

//!\brief The log of a difference of two log terms. (log_x - log_y)
// expm1 is more accurate than using exp if the difference is close to 0.
// std::log1p(-std::exp(difference)) = std::log(1 + (-std::exp(difference)))
//                                   = std::log(1 - std::exp(difference))
// std::log(-std::expm1(difference)) = std::log(-(std::exp(difference) - 1))
//                                   = std::log(1 - std::exp(difference))
[[nodiscard]] inline double substract(double const log_x, double const log_y) noexcept
{
    double const difference{log_y - log_x};
    return log_x + difference > -ln_2 ? std::log(-std::expm1(difference)) : std::log1p(-std::exp(difference));
}

struct add_fn
{
    [[nodiscard]] double operator()(double const log_x, double const log_y) const noexcept
    {
        return add(log_x, log_y);
    };
};

} // namespace raptor::logspace

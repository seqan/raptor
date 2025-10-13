// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides functionality for working with probabilities in logspace.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <algorithm> // for max
#include <cmath>     // for exp, log1p, expm1, log, abs
#include <cstdlib>   // for abs
#include <limits>    // for numeric_limits

namespace raptor::logspace
{

constexpr double ln_2{0.693147180559945};
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
    return add(add(log_x, log_y), logs...);
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

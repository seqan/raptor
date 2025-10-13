// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::threshold::pascal_row.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <cmath>   // for log
#include <cstddef> // for size_t
#include <vector>  // for vector

#include <raptor/threshold/pascal_row.hpp> // for pascal_row

namespace raptor::threshold
{

[[nodiscard]] std::vector<double> pascal_row(size_t const n)
{
    std::vector<double> result(n + 1);

    for (size_t i = 1; i <= n; ++i)
        result[i] = result[i - 1] + std::log((n + 1 - i) / i);

    return result;
}

} // namespace raptor::threshold

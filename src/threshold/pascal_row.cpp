// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Implements raptor::threshold::pascal_row.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <cmath>
#include <cstddef>

#include <raptor/threshold/pascal_row.hpp>

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

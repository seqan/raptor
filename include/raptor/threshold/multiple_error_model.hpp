// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::threshold::multiple_error_model.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <cstddef>
#include <vector>

namespace raptor::threshold
{

[[nodiscard]] std::vector<double> multiple_error_model(size_t const number_of_minimisers,
                                                       size_t const errors,
                                                       std::vector<double> const & affected_by_one_error_prob);

} // namespace raptor::threshold

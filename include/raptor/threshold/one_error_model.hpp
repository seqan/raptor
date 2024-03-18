// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::threshold::one_error_model.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <cstddef>
#include <vector>

namespace raptor::threshold
{

[[nodiscard]] std::vector<double> one_error_model(size_t const kmer_size,
                                                  double const p_mean,
                                                  std::vector<double> const & affected_by_one_error_indirectly_prob);

} // namespace raptor::threshold

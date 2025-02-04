// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::threshold::one_indirect_error_model.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <seqan3/search/kmer_index/shape.hpp>

namespace raptor::threshold
{

[[nodiscard]] std::vector<double>
one_indirect_error_model(size_t const query_length, size_t const window_size, seqan3::shape const shape);

} // namespace raptor::threshold

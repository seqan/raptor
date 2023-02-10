// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

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

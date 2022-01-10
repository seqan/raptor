// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <utility>
#include <vector>

namespace raptor::detail
{

std::tuple<double, std::vector<double>> simple_model(size_t const kmer_size,
                                                     std::vector<double> const & proba_x,
                                                     std::vector<double> const & indirect_errors);

} // namespace raptor::detail

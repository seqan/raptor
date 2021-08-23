// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <raptor/search/detail/helper.hpp>
#include <raptor/shared.hpp>

namespace raptor
{

void do_cerealisation_out(std::vector<size_t> const & vec, search_arguments const & arguments);
bool do_cerealisation_in(std::vector<size_t> & vec, search_arguments const & arguments);
std::vector<size_t> precompute_threshold(size_t const pattern_size,
                                         size_t const window_size,
                                         seqan3::shape const shape,
                                         size_t const errors,
                                         double const tau);

} // namespace raptor

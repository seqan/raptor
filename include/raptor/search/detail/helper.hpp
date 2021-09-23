// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <seqan3/search/kmer_index/shape.hpp>

namespace raptor::detail
{

std::vector<size_t> pascal_row(size_t n);

std::tuple<double, std::vector<double>>
simple_model(size_t const kmer_size, std::vector<double> const & proba_x, std::vector<double> const & indirect_errors);

void impl(size_t const minimizers_left,
          std::vector<double> const & proba,
          std::vector<size_t> error_distribution,
          size_t const current_error_index,
          double & result);

double enumerate_all_errors(size_t const number_of_minimizers, size_t const errors, std::vector<double> const & proba);

std::vector<double> destroyed_indirectly_by_error(size_t const pattern_size,
                                                  size_t const window_size,
                                                  seqan3::shape const shape);

} // namespace raptor::detail

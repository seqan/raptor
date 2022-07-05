// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <cmath>

#include <raptor/build/hibf/bin_size_in_bits.hpp>

namespace raptor::hibf
{

size_t bin_size_in_bits(build_arguments const & arguments, size_t const number_of_kmers_to_be_stored)
{
    double const numerator{-static_cast<double>(number_of_kmers_to_be_stored * arguments.hash)};
    double const denominator{std::log(1 - std::exp(std::log(arguments.fpr) / arguments.hash))};
    double const result{std::ceil(numerator / denominator)};
    return result;
}

} // namespace raptor::hibf

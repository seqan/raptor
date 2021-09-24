// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <cmath>

#include <raptor/build/hibf/bin_size_in_bits.hpp>

namespace raptor::hibf
{

size_t bin_size_in_bits(build_config const & config, size_t const number_of_kmers_to_be_stored)
{
    return std::ceil( - static_cast<double>(number_of_kmers_to_be_stored * config.hash_funs) /
                     std::log(1 - std::exp(std::log(config.FPR) / config.hash_funs)));
}

} // namespace raptor::hibf

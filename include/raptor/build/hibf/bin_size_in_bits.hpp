// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <raptor/argument_parsing/build_arguments.hpp>

namespace raptor::hibf
{

size_t bin_size_in_bits(build_arguments const & arguments, size_t const number_of_kmers_to_be_stored);

} // namespace raptor::hibf

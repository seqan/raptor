// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <robin_hood.h>

#include <raptor/argument_parsing/build_arguments.hpp>
#include <raptor/build/hibf/build_data.hpp>

namespace raptor::hibf
{

template <seqan3::data_layout data_layout_mode>
size_t initialise_max_bin_kmers(robin_hood::unordered_flat_set<size_t> & kmers,
                                std::vector<int64_t> & ibf_positions,
                                std::vector<int64_t> & filename_indices,
                                lemon::ListDigraph::Node const & node,
                                build_data<data_layout_mode> & data,
                                build_arguments const & arguments);

} // namespace raptor::hibf

// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides raptor::hibf::hierarchical_build.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <robin_hood.h>

#include <raptor/build/hibf/build_data.hpp>

namespace raptor::hibf
{

size_t hierarchical_build(lemon::ListDigraph::Node const & root_node, build_data & data);

size_t hierarchical_build(robin_hood::unordered_flat_set<size_t> & parent_kmers,
                          lemon::ListDigraph::Node const & current_node,
                          build_data & data,
                          bool is_root);

} // namespace raptor::hibf

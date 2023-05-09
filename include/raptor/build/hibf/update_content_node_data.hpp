// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides raptor::hibf::update_content_node_data.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <chopper/layout/layout.hpp>

#include <raptor/build/hibf/node_data.hpp>

namespace raptor::hibf
{

void update_content_node_data(std::vector<chopper::layout::layout::user_bin> && layout_user_bins,
                              lemon::ListDigraph & ibf_graph,
                              lemon::ListDigraph::NodeMap<node_data> & node_map);

} // namespace raptor::hibf

// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Implements raptor::hibf::update_header_node_data.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#include <algorithm>
#include <cassert>

#include <raptor/build/hibf/update_header_node_data.hpp>

namespace raptor::hibf
{

void update_header_node_data(std::vector<chopper::layout::layout::max_bin> & header_max_bins,
                             lemon::ListDigraph & ibf_graph,
                             lemon::ListDigraph::NodeMap<node_data> & node_map)
{
    // sort records ascending by the number of bin indices (corresponds to the IBF levels)
    std::ranges::sort(header_max_bins,
                      [](auto const & r, auto const & l)
                      {
                          return r.previous_TB_indices.size() < l.previous_TB_indices.size();
                      });

    for (auto const & [bin_indices, max_id] : header_max_bins)
    {
        // we assume that the header lines are in the correct order
        // go down the tree until you find the matching parent
        auto it = bin_indices.begin();
        lemon::ListDigraph::Node current_node = ibf_graph.nodeFromId(0); // start at root

        while (it != (bin_indices.end() - 1))
        {
            for (lemon::ListDigraph::OutArcIt arc_it(ibf_graph, current_node); arc_it != lemon::INVALID; ++arc_it)
            {
                auto target = ibf_graph.target(arc_it);
                if (node_map[target].parent_bin_index == *it)
                {
                    current_node = target;
                    break;
                }
            }
            ++it;
        }

        auto new_node = ibf_graph.addNode();
        ibf_graph.addArc(current_node, new_node);
        node_map.set(new_node, {bin_indices.back(), max_id, 0, lemon::INVALID, {}});

        if (node_map[current_node].max_bin_index == bin_indices.back())
            node_map[current_node].favourite_child = new_node;
    }
}

} // namespace raptor::hibf

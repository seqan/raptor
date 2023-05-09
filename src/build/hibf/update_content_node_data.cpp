// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Implements raptor::hibf::update_content_node_data.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#include <cassert>

#include <raptor/build/hibf/update_content_node_data.hpp>

namespace raptor::hibf
{

void update_content_node_data(std::vector<chopper::layout::layout::user_bin> & layout_user_bins,
                              lemon::ListDigraph & ibf_graph,
                              lemon::ListDigraph::NodeMap<node_data> & node_map)
{
    // parse lines
    // -------------------------------------------------------------------------
    for (size_t user_bin = 0; user_bin < layout_user_bins.size(); ++user_bin)
    {
        auto const & record = layout_user_bins[user_bin];

        // go down the tree until you find the matching parent
        lemon::ListDigraph::Node current_node = ibf_graph.nodeFromId(0); // start at root

        for (size_t i = 0; i < record.previous_TB_indices.size(); ++i)
        {
            size_t const bin = record.previous_TB_indices[i];
            size_t const num_tbs = 1;
            auto & current_data = node_map[current_node];

            // update number of technical bins in current_node-IBF
            current_data.number_of_technical_bins = std::max(current_data.number_of_technical_bins, bin + num_tbs);

#ifndef NDEBUG
            bool found_next_node{false}; // sanity check
#endif
            for (lemon::ListDigraph::OutArcIt arc_it(ibf_graph, current_node); arc_it != lemon::INVALID; ++arc_it)
            {
                auto target = ibf_graph.target(arc_it);
                if (node_map[target].parent_bin_index == bin)
                {
                    current_node = target;
#ifndef NDEBUG
                    found_next_node = true;
#endif
                    break;
                }
            }
            assert(found_next_node);
        }

        size_t const bin = record.storage_TB_id;
        size_t const num_tbs = record.number_of_technical_bins;
        auto & current_data = node_map[current_node];

        // update number of technical bins in current_node-IBF
        current_data.number_of_technical_bins = std::max(current_data.number_of_technical_bins, bin + num_tbs);

        if (record.storage_TB_id == current_data.max_bin_index)
            current_data.remaining_records.insert(current_data.remaining_records.begin(), record);
        else
            current_data.remaining_records.push_back(record);
    }
}

} // namespace raptor::hibf

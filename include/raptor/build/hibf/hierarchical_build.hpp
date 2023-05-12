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

#include <lemon/list_graph.h> /// Must be first include.

#include <robin_hood.h>

#include <raptor/build/hibf/build_data.hpp>
#include <raptor/build/hibf/construct_ibf.hpp>
#include <raptor/build/hibf/hierarchical_build.hpp>
#include <raptor/build/hibf/insert_into_ibf.hpp>
#include <raptor/build/hibf/loop_over_children.hpp>
#include <raptor/build/hibf/update_parent_kmers.hpp>
#include <raptor/build/hibf/update_user_bins.hpp>

namespace raptor::hibf
{

template <typename input_range_type>
size_t hierarchical_build(robin_hood::unordered_flat_set<size_t> & parent_kmers,
                          lemon::ListDigraph::Node const & current_node,
                          build_data<input_range_type> & data,
                          bool is_root)
{
    auto & current_node_data = data.node_map[current_node];

    size_t const ibf_pos{data.request_ibf_idx()};

    std::vector<int64_t> ibf_positions(current_node_data.number_of_technical_bins, ibf_pos);
    std::vector<int64_t> filename_indices(current_node_data.number_of_technical_bins, -1);

    // initialize lower level IBF
    size_t max_bin_tbs;
    seqan3::interleaved_bloom_filter ibf;

    if (auto & node_data = data.node_map[node]; node_data.favourite_child != lemon::INVALID) // max bin is a merged bin
    {
        // recursively initialize favourite child first
        robin_hood::unordered_flat_set<size_t> kmers{};
        ibf_positions[node_data.max_bin_index] = hierarchical_build(kmers, node_data.favourite_child, data, false);

        ibf = construct_ibf(parent_kmers, kmers, max_bin_tbs, current_node, data, is_root);
        max_bin_tbs = 1;
    }
    else // max bin is not a merged bin
    {
        // we assume that the max record is at the beginning of the list of remaining records.
        auto const & record = node_data.remaining_records[0];
        update_user_bins(filename_indices, record);

        ibf = construct_ibf(parent_kmers, data.input[record.idx], max_bin_tbs, current_node, data, is_root);
        max_bin_tbs = record.number_of_technical_bins;
    }

    // parse all other children (merged bins) of the current ibf
    loop_over_children(parent_kmers, ibf, ibf_positions, current_node, data, is_root);

    // If max bin was a merged bin, process all remaining records, otherwise the first one has already been processed
    size_t const start{(current_node_data.favourite_child != lemon::INVALID) ? 0u : 1u};
    for (size_t i = start; i < current_node_data.remaining_records.size(); ++i)
    {
        auto const & record = current_node_data.remaining_records[i];

        if (is_root && record.number_of_technical_bins == 1) // no splitting needed
        {
            insert_into_ibf(data, record, ibf);
        }
        else
        {
            insert_into_ibf(data.input[record.idx],
                            record.number_of_technical_bins,
                            record.storage_TB_id,
                            ibf,
                            data.arguments.fill_ibf_timer);
            if (!is_root)
                update_parent_kmers(parent_kmers, kmers, data.arguments.merge_kmers_timer);
        }

        update_user_bins(filename_indices, record);
        kmers.clear();
    }

    data.hibf.ibf_vector[ibf_pos] = std::move(ibf);
    data.hibf.next_ibf_id[ibf_pos] = std::move(ibf_positions);
    data.hibf.user_bins.bin_indices_of_ibf(ibf_pos) = std::move(filename_indices);

    return ibf_pos;
}

} // namespace raptor::hibf

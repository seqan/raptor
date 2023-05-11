// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Implements raptor::hibf::construct_ibf.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <lemon/list_graph.h> /// Must be first include.

#include <raptor/build/hibf/bin_size_in_bits.hpp>
#include <raptor/build/hibf/construct_ibf.hpp>
#include <raptor/build/hibf/insert_into_ibf.hpp>
#include <raptor/build/hibf/update_parent_kmers.hpp>

namespace raptor::hibf
{

seqan3::interleaved_bloom_filter<> construct_ibf(robin_hood::unordered_flat_set<size_t> & parent_kmers,
                                                 robin_hood::unordered_flat_set<size_t> & kmers,
                                                 size_t const number_of_bins,
                                                 lemon::ListDigraph::Node const & node,
                                                 build_data & data,
                                                 bool is_root)
{
    auto & node_data = data.node_map[node];

    size_t const kmers_per_bin{static_cast<size_t>(std::ceil(static_cast<double>(kmers.size()) / number_of_bins))};
    double const bin_bits{static_cast<double>(bin_size_in_bits(data.arguments, kmers_per_bin))};
    seqan3::bin_size const bin_size{static_cast<size_t>(std::ceil(bin_bits * data.fp_correction[number_of_bins]))};
    seqan3::bin_count const bin_count{node_data.number_of_technical_bins};

    timer<concurrent::no> local_index_allocation_timer{};
    local_index_allocation_timer.start();
    seqan3::interleaved_bloom_filter<> ibf{bin_count, bin_size, seqan3::hash_function_count{data.arguments.hash}};
    local_index_allocation_timer.stop();
    data.arguments.index_allocation_timer += local_index_allocation_timer;

    insert_into_ibf(kmers, number_of_bins, node_data.max_bin_index, ibf, data.arguments.fill_ibf_timer);
    if (!is_root)
        update_parent_kmers(parent_kmers, kmers, data.arguments.merge_kmers_timer);

    return ibf;
}

} // namespace raptor::hibf

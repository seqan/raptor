// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <lemon/list_graph.h> /// Must be first include.

#include <raptor/build/hibf/bin_size_in_bits.hpp>
#include <raptor/build/hibf/construct_ibf.hpp>
#include <raptor/build/hibf/insert_into_ibf.hpp>

namespace raptor::hibf
{

template <seqan3::data_layout data_layout_mode>
seqan3::interleaved_bloom_filter<> construct_ibf(robin_hood::unordered_flat_set<size_t> & parent_kmers,
                                                 robin_hood::unordered_flat_set<size_t> & kmers,
                                                 size_t const number_of_bins,
                                                 lemon::ListDigraph::Node const & node,
                                                 build_data<data_layout_mode> & data,
                                                 build_arguments const & arguments,
                                                 bool is_root)
{
    auto & node_data = data.node_map[node];

    size_t const kmers_per_bin{static_cast<size_t>(std::ceil(static_cast<double>(kmers.size()) / number_of_bins))};
    double const bin_bits{static_cast<double>(bin_size_in_bits(arguments, kmers_per_bin))};
    seqan3::bin_size const bin_size{static_cast<size_t>(std::ceil(bin_bits * data.fp_correction[number_of_bins]))};
    seqan3::bin_count const bin_count{node_data.number_of_technical_bins};
    seqan3::interleaved_bloom_filter<> ibf{bin_count, bin_size, seqan3::hash_function_count{arguments.hash}};

    insert_into_ibf(parent_kmers, kmers, number_of_bins, node_data.max_bin_index, ibf, is_root);

    return ibf;
}

template seqan3::interleaved_bloom_filter<>
construct_ibf<seqan3::data_layout::uncompressed>(robin_hood::unordered_flat_set<size_t> &,
                                                 robin_hood::unordered_flat_set<size_t> &,
                                                 size_t const,
                                                 lemon::ListDigraph::Node const &,
                                                 build_data<seqan3::data_layout::uncompressed> &,
                                                 build_arguments const &,
                                                 bool);

template seqan3::interleaved_bloom_filter<>
construct_ibf<seqan3::data_layout::compressed>(robin_hood::unordered_flat_set<size_t> &,
                                               robin_hood::unordered_flat_set<size_t> &,
                                               size_t const,
                                               lemon::ListDigraph::Node const &,
                                               build_data<seqan3::data_layout::compressed> &,
                                               build_arguments const &,
                                               bool);

} // namespace raptor::hibf

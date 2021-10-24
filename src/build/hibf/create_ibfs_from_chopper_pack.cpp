// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <lemon/list_graph.h> /// Must be first include.

#include <raptor/build/hibf/create_ibfs_from_chopper_pack.hpp>
#include <raptor/build/hibf/hierarchical_build.hpp>
#include <raptor/build/hibf/read_chopper_pack_file.hpp>

namespace raptor::hibf
{

void create_ibfs_from_chopper_pack(build_data & data, build_arguments const & arguments)
{
    // std::cerr << "[DEBUG] create_ibfs_from_chopper_pack\n"
    //           << "        &data:          " << &data << '\n'
    //           << "        &data.node_map: " << &data.node_map << "\n\n";
    read_chopper_pack_file(data, arguments.bin_file);
    lemon::ListDigraph::Node root = data.ibf_graph.nodeFromId(0); // root node = high level IBF node
    robin_hood::unordered_flat_set<size_t> root_kmers{};

    size_t const t_max{data.node_map[root].number_of_technical_bins};
    data.compute_fp_correction(t_max, arguments.hash, arguments.fpr);

    hierarchical_build(root_kmers, root, data, arguments, true);
}

} // namespace raptor::hibf

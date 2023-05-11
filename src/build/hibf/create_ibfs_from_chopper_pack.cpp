// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Implements raptor::hibf::create_ibfs_from_chopper_pack.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <lemon/list_graph.h> /// Must be first include.

#include <chopper/layout/compute_fp_correction.hpp>

#include <raptor/build/hibf/create_ibfs_from_chopper_pack.hpp>
#include <raptor/build/hibf/hierarchical_build.hpp>
#include <raptor/build/hibf/initialise_build_tree.hpp>
#include <raptor/build/hibf/read_chopper_pack_file.hpp>

namespace raptor::hibf
{

void create_ibfs_from_chopper_pack(build_data & data, build_arguments const & arguments)
{
    chopper::layout::layout hibf_layout = read_chopper_pack_file(data.filenames, arguments.bin_file);

    size_t const number_of_ibfs = hibf_layout.max_bins.size() + 1;

    data.hibf.ibf_vector.resize(number_of_ibfs);
    data.hibf.user_bins.set_ibf_count(number_of_ibfs);
    data.hibf.user_bins.set_user_bin_count(hibf_layout.user_bins.size());
    data.hibf.next_ibf_id.resize(number_of_ibfs);

    initialise_build_tree(hibf_layout, data.ibf_graph, data.node_map);

    lemon::ListDigraph::Node root = data.ibf_graph.nodeFromId(0); // root node = high level IBF node
    robin_hood::unordered_flat_set<size_t> root_kmers{};

    size_t const t_max{data.node_map[root].number_of_technical_bins};
    data.fp_correction = chopper::layout::compute_fp_correction(arguments.fpr, arguments.hash, t_max);

    hierarchical_build(root_kmers, root, data, arguments, true);
}

} // namespace raptor::hibf

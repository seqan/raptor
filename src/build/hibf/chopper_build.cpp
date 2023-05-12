// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Implements raptor::hibf::chopper_build.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <lemon/list_graph.h> /// Must be first include.

#include <chopper/layout/compute_fp_correction.hpp>

#include <raptor/build/hibf/build_data.hpp>
#include <raptor/build/hibf/chopper_build.hpp>
#include <raptor/build/hibf/hierarchical_build.hpp>
#include <raptor/build/hibf/initialise_build_tree.hpp>
#include <raptor/build/hibf/input_base.hpp>
#include <raptor/build/hibf/read_chopper_pack_file.hpp>
#include <raptor/build/store_index.hpp>
#include <raptor/file_reader.hpp>
#include <raptor/index.hpp>

namespace raptor
{

class input_derived : public raptor::hibf::input_base
{
public:
    void hash_into(size_t const user_bin_number, robin_hood::unordered_flat_set<size_t> & to) const override
    {
        if (arguments.input_is_minimiser)
        {
            file_reader<file_types::minimiser> const reader{};
            reader.hash_into(filenames[user_bin_number], std::inserter(to, to.begin()));
        }
        else
        {
            file_reader<file_types::sequence> const reader{arguments.shape, arguments.window_size};
            reader.hash_into(filenames[user_bin_number], std::inserter(to, to.begin()));
        }
    }

    build_arguments const & arguments{};

    std::vector<std::vector<std::string>> const & filenames{};

    input_derived(build_arguments const & args, std::vector<std::vector<std::string>> const & filenames_) :
        arguments{args},
        filenames{filenames_}
    {}
};

} // namespace raptor

namespace raptor::hibf
{

void chopper_build(build_arguments const & arguments)
{
    std::vector<std::vector<std::string>> filenames{};

    chopper::layout::layout hibf_layout = read_chopper_pack_file(filenames, arguments.bin_file);

    build_data data{.arguments = arguments, .input = input_derived{arguments, filenames}};

    size_t const number_of_ibfs = hibf_layout.max_bins.size() + 1;
    data.hibf.ibf_vector.resize(number_of_ibfs);
    data.hibf.user_bins.set_ibf_count(number_of_ibfs);
    data.hibf.user_bins.set_user_bin_count(hibf_layout.user_bins.size());
    data.hibf.next_ibf_id.resize(number_of_ibfs);

    initialise_build_tree(hibf_layout, data.ibf_graph, data.node_map);
    lemon::ListDigraph::Node root_node = data.ibf_graph.nodeFromId(0); // root node = top-level IBF node

    size_t const t_max{data.node_map[root_node].number_of_technical_bins};
    data.fp_correction = chopper::layout::compute_fp_correction(arguments.fpr, arguments.hash, t_max);

    hierarchical_build(root_node, data);

    arguments.index_allocation_timer.start();
    raptor_index<hierarchical_interleaved_bloom_filter> index{window{arguments.window_size},
                                                              arguments.shape,
                                                              arguments.parts,
                                                              filenames,
                                                              arguments.fpr,
                                                              std::move(data.hibf)};
    arguments.index_allocation_timer.stop();

    arguments.store_index_timer.start();
    store_index(arguments.out_path, std::move(index), arguments);
    arguments.store_index_timer.stop();
}

} // namespace raptor::hibf

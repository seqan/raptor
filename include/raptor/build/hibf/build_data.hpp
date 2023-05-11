// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides raptor::hibf::build_data.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <lemon/list_graph.h> /// Must be first include.

#include <atomic>

#include <raptor/build/hibf/node_data.hpp>
#include <raptor/hierarchical_interleaved_bloom_filter.hpp>

namespace raptor::hibf
{

struct build_data
{
    std::atomic<size_t> ibf_number{};

    size_t number_of_user_bins{};
    size_t number_of_ibfs{};

    std::vector<std::vector<std::string>> filenames{};

    lemon::ListDigraph ibf_graph{};
    lemon::ListDigraph::NodeMap<node_data> node_map{ibf_graph};

    hierarchical_interleaved_bloom_filter hibf{};
    std::vector<double> fp_correction{};

    size_t request_ibf_idx()
    {
        return std::atomic_fetch_add(&ibf_number, 1u);
    }

    void resize()
    {
        hibf.ibf_vector.resize(number_of_ibfs);
        hibf.user_bins.set_ibf_count(number_of_ibfs);
        hibf.user_bins.set_user_bin_count(number_of_user_bins);
        hibf.next_ibf_id.resize(number_of_ibfs);
    }
};

} // namespace raptor::hibf

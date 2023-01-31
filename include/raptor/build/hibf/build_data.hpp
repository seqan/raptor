// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <atomic>
#include <seqan3/std/new>

#include <raptor/build/hibf/node_data.hpp>
#include <raptor/hierarchical_interleaved_bloom_filter.hpp>

namespace raptor::hibf
{

template <seqan3::data_layout data_layout_mode>
struct build_data
{
    alignas(std::hardware_destructive_interference_size) std::atomic<size_t> ibf_number{};
    alignas(std::hardware_destructive_interference_size) std::atomic<size_t> user_bin_number{};

    size_t number_of_user_bins{};
    size_t number_of_ibfs{};

    lemon::ListDigraph ibf_graph{};
    lemon::ListDigraph::NodeMap<node_data> node_map{ibf_graph};

    hierarchical_interleaved_bloom_filter<data_layout_mode> hibf{};
    std::vector<double> fp_correction{};

    size_t request_ibf_idx()
    {
        return std::atomic_fetch_add(&ibf_number, 1u);
    }

    size_t request_user_bin_idx()
    {
        return std::atomic_fetch_add(&user_bin_number, 1u);
    }

    void resize()
    {
        hibf.ibf_vector.resize(number_of_ibfs);
        hibf.user_bins.set_ibf_count(number_of_ibfs);
        hibf.user_bins.set_user_bin_count(number_of_user_bins);
        hibf.next_ibf_id.resize(number_of_ibfs);
    }

    void compute_fp_correction(size_t const tmax, size_t const hash, double const fpr)
    {
        fp_correction.resize(tmax + 1, 1.0);

        double const denominator = std::log(1 - std::exp(std::log(fpr) / hash));

        for (size_t i = 2; i <= tmax; ++i)
        {
            double const tmp = 1.0 - std::pow(1 - fpr, static_cast<double>(i));
            fp_correction[i] = std::log(1 - std::exp(std::log(tmp) / hash)) / denominator;
            assert(fp_correction[i] >= 1.0);
        }
    }
};

} // namespace raptor::hibf

// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <seqan3/std/filesystem>
#include <vector>

#include <seqan3/search/kmer_index/shape.hpp>

#include <raptor/strong_types.hpp>

namespace raptor
{

struct build_arguments
{
    // Related to k-mers
    uint8_t kmer_size{20u};
    uint32_t window_size{kmer_size};
    window window_size_strong{kmer_size};
    std::string shape_string{};
    seqan3::shape shape{seqan3::ungapped{kmer_size}};
    bool compute_minimiser{false};
    bool disable_cutoffs{false};

    // Related to IBF
    std::filesystem::path out_path{"./"};
    std::string size{"1k"};
    uint64_t bins{64};
    uint64_t bits{4096};
    uint64_t hash{2};
    uint8_t parts{1u};
    double fpr{0.05};
    bool compressed{false};

    // General arguments
    std::vector<std::vector<std::string>> bin_path{};
    std::filesystem::path bin_file{};
    uint8_t threads{1u};
    bool is_socks{false};
    bool is_hibf{false};
};

} // namespace raptor

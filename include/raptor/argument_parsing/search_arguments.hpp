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

namespace raptor
{

struct search_arguments
{
    // Related to k-mers
    uint32_t window_size{20u};
    seqan3::shape shape{seqan3::ungapped{20u}};
    uint8_t shape_size{shape.size()};
    uint8_t shape_weight{shape.count()};
    uint8_t threads{1u};
    uint8_t parts{1u};

    // Related to thresholding
    double tau{0.99};
    double threshold{};
    uint64_t pattern_size{};
    uint8_t errors{0};
    bool treshold_was_set{false};

    // Related to IBF
    std::filesystem::path index_file{};
    bool compressed{false};

    // General arguments
    std::vector<std::vector<std::string>> bin_path{};
    std::filesystem::path query_file{};
    std::filesystem::path out_file{"search.out"};
    bool write_time{false};
    bool is_socks{false};
    bool is_hibf{false};
};

} // namespace raptor

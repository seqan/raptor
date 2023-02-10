// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <filesystem>
#include <vector>

#include <seqan3/search/kmer_index/shape.hpp>

namespace raptor
{

struct upgrade_arguments
{
    uint32_t window_size{};
    seqan3::shape shape{};
    bool compressed{};
    bool input_is_minimiser{};
    uint8_t parts{1u};
    uint8_t threads{1u};
    double fpr{std::numeric_limits<double>::quiet_NaN()};

    std::filesystem::path bin_file{};
    std::filesystem::path index_file{};
    std::filesystem::path output_file{};

    std::vector<std::vector<std::string>> bin_path{};
};

} // namespace raptor

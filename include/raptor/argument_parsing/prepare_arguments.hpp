// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <filesystem>
#include <vector>

#include <seqan3/search/kmer_index/shape.hpp>

#include <raptor/strong_types.hpp>

namespace raptor
{

struct prepare_arguments
{
    uint8_t kmer_size{20u};
    uint32_t window_size{kmer_size};
    std::string shape_string{};
    seqan3::shape shape{seqan3::ungapped{kmer_size}};
    bool use_filesize_dependent_cutoff{false};
    uint8_t kmer_count_cutoff{1u};

    std::filesystem::path out_dir{"./"};

    std::vector<std::vector<std::string>> bin_path{};
    std::filesystem::path bin_file{};
    uint8_t threads{1u};
};

} // namespace raptor

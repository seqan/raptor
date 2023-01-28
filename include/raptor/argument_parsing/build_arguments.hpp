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

#include <raptor/argument_parsing/prepare_arguments.hpp>
#include <raptor/strong_types.hpp>

namespace raptor
{

struct build_arguments
{
    // Related to k-mers
    uint8_t kmer_size{20u};
    uint32_t window_size{kmer_size};
    std::string shape_string{};
    seqan3::shape shape{seqan3::ungapped{kmer_size}};

    // Related to IBF
    std::filesystem::path out_path{"./"};
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
    bool is_hibf{false};
    bool input_is_minimiser{false};

    prepare_arguments make_prepare_arguments(std::filesystem::path const out_dir) const noexcept
    {
        return {.kmer_size{kmer_size},
                .window_size{window_size},
                .shape{shape},
                .use_filesize_dependent_cutoff{false},
                .kmer_count_cutoff{1u},
                .out_dir{out_dir},
                .bin_path{bin_path},
                .threads{threads}};
    }
};

} // namespace raptor

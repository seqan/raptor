// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <seqan3/std/filesystem>
#include <vector>

#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/kmer_index/shape.hpp>

#include <raptor/adjust_seed.hpp>

namespace raptor
{

//!\brief Strong type for passing the window size.
struct window { uint64_t v; };
//!\brief Strong type for passing the kmer size.
struct kmer { uint8_t v; };
//!\brief Strong type for passing number of bins.
struct bins { uint64_t v; };
//!\brief Strong type for passing number of bits.
struct bits { uint64_t v; };
//!\brief Strong type for passing number of hash functions.
struct hashes { uint64_t v; };

struct dna4_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

struct build_arguments
{
    // Related to k-mers
    uint32_t window_size{20u};
    uint8_t kmer_size{20u};
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
    bool compressed{false};

    // General arguments
    std::vector<std::vector<std::string>> bin_path{};
    std::filesystem::path bin_file{};
    uint8_t threads{1u};
    bool is_socks{false};
};

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
};

struct upgrade_arguments
{
    uint32_t window_size{};
    uint8_t kmer_size{};
    seqan3::shape shape{};
    uint8_t parts{1u};
    bool compressed{false};

    std::filesystem::path bin_file{};
    std::filesystem::path in_file{};
    std::filesystem::path out_file{};

    std::vector<std::vector<std::string>> bin_path{};
};

} // namespace raptor

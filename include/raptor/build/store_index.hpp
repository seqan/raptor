// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <seqan3/std/filesystem>

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

#include <raptor/index.hpp>
#include <raptor/strong_types.hpp>

namespace raptor
{

template <typename data_t, typename arguments_t>
static inline void store_index(std::filesystem::path const & path,
                               raptor_index<data_t> const & index,
                               arguments_t const & arguments)
{
    std::ofstream os{path, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(index);
}

template <seqan3::data_layout layout, typename arguments_t>
static inline void store_index(std::filesystem::path const & path,
                               seqan3::interleaved_bloom_filter<layout> && ibf,
                               arguments_t const & arguments)
{
    raptor_index<seqan3::interleaved_bloom_filter<layout>> index{window{arguments.window_size},
                                                                 arguments.shape,
                                                                 arguments.parts,
                                                                 arguments.compressed,
                                                                 arguments.bin_path,
                                                                 std::move(ibf)};

    std::ofstream os{path, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(index);
}

} // namespace raptor


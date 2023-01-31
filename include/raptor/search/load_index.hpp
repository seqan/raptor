// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <chrono>

#include <raptor/argument_parsing/search_arguments.hpp>
#include <raptor/index.hpp>

namespace raptor
{

template <typename index_t>
void load_index(index_t & index, search_arguments const & arguments, size_t const part, double & index_io_time)
{
    std::filesystem::path index_file{arguments.index_file};
    index_file += "_" + std::to_string(part);

    load_index(index, index_file, index_io_time);
}

template <typename index_t>
void load_index(index_t & index, search_arguments const & arguments, double & index_io_time)
{
    load_index(index, arguments.index_file, index_io_time);
}

template <typename index_t>
void load_index(index_t & index, std::filesystem::path const & path, double & index_io_time)
{
    std::ifstream is{path, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};

    auto start = std::chrono::high_resolution_clock::now();
    iarchive(index);
    auto end = std::chrono::high_resolution_clock::now();

    index_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
}

} // namespace raptor

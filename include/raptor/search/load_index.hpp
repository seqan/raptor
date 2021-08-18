// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <chrono>
#include <seqan3/std/filesystem>

#include <raptor/index.hpp>
#include <raptor/shared.hpp>

namespace raptor
{

template <typename t>
void load_index(t & index, search_arguments const & arguments, size_t const part, double & index_io_time)
{
    std::filesystem::path index_file{arguments.index_file};
    index_file += "_" + std::to_string(part);

    std::ifstream is{index_file, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};

    auto start = std::chrono::high_resolution_clock::now();
    iarchive(index);
    auto end = std::chrono::high_resolution_clock::now();

    index_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
}

template <typename t>
void load_index(t & index, search_arguments const & arguments, double & index_io_time)
{
    std::ifstream is{arguments.index_file, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};

    auto start = std::chrono::high_resolution_clock::now();
    iarchive(index);
    auto end = std::chrono::high_resolution_clock::now();

    index_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
}

} // namespace raptor

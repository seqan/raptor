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

namespace detail
{

template <typename index_t>
void load_index(index_t & index, std::filesystem::path const & path)
{
    std::ifstream is{path, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};

    iarchive(index);
}

} // namespace detail

template <typename index_t>
void load_index(index_t & index, search_arguments const & arguments, size_t const part)
{
    std::filesystem::path index_file{arguments.index_file};
    index_file += "_" + std::to_string(part);
    arguments.load_index_timer.start();
    detail::load_index(index, index_file);
    arguments.load_index_timer.stop();
}

template <typename index_t>
void load_index(index_t & index, search_arguments const & arguments)
{
    arguments.load_index_timer.start();
    detail::load_index(index, arguments.index_file);
    arguments.load_index_timer.stop();
}

} // namespace raptor

// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::load_index.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <chrono>
#include <fstream>

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

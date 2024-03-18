// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::store_index.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <filesystem>
#include <fstream>

#include <raptor/index.hpp>
#include <raptor/strong_types.hpp>

namespace raptor
{

template <typename data_t>
static inline void store_index(std::filesystem::path const & path, raptor_index<data_t> && index)
{
    std::ofstream os{path, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(index);
}

} // namespace raptor

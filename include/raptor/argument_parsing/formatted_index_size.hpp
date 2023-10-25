// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides raptor::formatted_index_size.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <cassert>
#include <cstddef>
#include <filesystem>
#include <string>

#include <raptor/argument_parsing/formatted_bytes.hpp>

namespace raptor
{

[[nodiscard]] inline size_t index_size_in_KiB(std::filesystem::path index_path, uint8_t const parts)
{
    size_t index_size_in_bytes{};
    if (parts == 1u)
    {
        index_size_in_bytes = std::filesystem::file_size(index_path);
    }
    else
    {
        for (size_t part = 0u; part < parts; ++part)
        {
            index_size_in_bytes += std::filesystem::file_size(index_path.string() + "_" + std::to_string(part));
        }
    }

    return index_size_in_bytes >> 10;
}

[[nodiscard]] inline std::string formatted_index_size(std::filesystem::path index_path, uint8_t const parts)
{
    return formatted_bytes(index_size_in_KiB(index_path, parts) << 10);
}

} // namespace raptor

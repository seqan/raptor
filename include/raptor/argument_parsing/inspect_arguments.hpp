// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::inspect_arguments.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <cstdint>
#include <filesystem>

namespace raptor
{

//!\brief The output formats of `raptor inspect`.
enum class inspect_format : uint8_t
{
    text, //!< `## Key = value` lines, like the header of a `raptor search` result.
    json  //!< A single JSON object.
};

struct inspect_arguments
{
    //!\brief The index to inspect. For a partitioned index, the path without the `_<part>` suffix.
    std::filesystem::path index_file{};
    //!\brief The file to write the metadata to. Empty for `std::cout`.
    std::filesystem::path out_file{};
    inspect_format format{inspect_format::text};
    //!\brief Whether to list the files of every user bin.
    bool user_bins{false};
    //!\brief Whether to load the whole index to also report its structure.
    bool detailed{false};

    // Determined by the parser, not set on the command line.
    //!\brief Whether `index_file` names a partitioned index, i.e. `<index_file>_0` exists.
    bool is_partitioned{false};
};

} // namespace raptor

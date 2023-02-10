// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides layout file string prefixes.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <string_view>

namespace raptor::hibf
{

constexpr std::string_view pack_file_first_line_prefix{"#HIGH_LEVEL_IBF"};
constexpr std::string_view hibf_prefix{"HIGH_LEVEL_IBF"};
static_assert(pack_file_first_line_prefix.ends_with(hibf_prefix));
constexpr std::string_view merged_bin_prefix{"MERGED_BIN"};
constexpr std::string_view split_bin_prefix{"SPLIT_BIN"};
constexpr size_t merged_bin_prefix_length{merged_bin_prefix.size()};

} // namespace raptor::hibf

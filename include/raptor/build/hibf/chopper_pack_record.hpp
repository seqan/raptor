// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides raptor::hibf::chopper_pack_record.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <string>
#include <tuple>
#include <vector>

#include <chopper/layout/layout.hpp>

namespace raptor::hibf
{

struct chopper_pack_record
{
    std::vector<std::string> filenames{};
    chopper::layout::layout::user_bin user_bin_info{};

    bool operator==(chopper_pack_record const & other) const
    {
        return std::tie(filenames, user_bin_info) == std::tie(other.filenames, other.user_bin_info);
    }

    bool operator!=(chopper_pack_record const & other) const
    {
        return std::tie(filenames, user_bin_info) != std::tie(other.filenames, other.user_bin_info);
    }
};

} // namespace raptor::hibf

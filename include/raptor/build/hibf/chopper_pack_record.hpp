// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <string>
#include <tuple>
#include <vector>

namespace raptor::hibf
{

struct chopper_pack_record
{
    std::vector<std::string> filenames{};
    std::vector<size_t> bin_indices{};
    std::vector<size_t> number_of_bins{};
    std::vector<size_t> estimated_sizes{};

    bool operator==(chopper_pack_record const & other) const
    {
        return std::tie(filenames, bin_indices, number_of_bins, estimated_sizes)
            == std::tie(other.filenames, other.bin_indices, other.number_of_bins, other.estimated_sizes);
    }

    bool operator!=(chopper_pack_record const & other) const
    {
        return std::tie(filenames, bin_indices, number_of_bins, estimated_sizes)
            != std::tie(other.filenames, other.bin_indices, other.number_of_bins, other.estimated_sizes);
    }
};

} // namespace raptor::hibf

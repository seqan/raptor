// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <bit>
#include <cstddef>

namespace raptor
{

struct partition_config
{
    partition_config(partition_config const &) = default;
    partition_config(partition_config &&) = default;
    partition_config & operator=(partition_config const &) = default;
    partition_config & operator=(partition_config &&) = default;
    ~partition_config() = default;

    explicit partition_config(size_t const parts) : partitions{parts}
    {
        size_t const next_power = next_power_of_four(partitions);
        suffixes_per_part = next_power / partitions;
        mask = next_power - 1;
    }

    size_t partitions{};
    size_t suffixes_per_part{};
    size_t mask{};

    static constexpr size_t next_power_of_four(size_t number)
    {
        if (number == 0ULL || number == 1ULL)
            return 1ULL;

        --number;
        int const highest_set_bit_pos = std::bit_width(number);
        int const shift_amount = (highest_set_bit_pos + (highest_set_bit_pos & 1)) - 2;
        //                       (           Next multiple of two                )   4 has two zeros
        return 0b0100ULL << shift_amount;
    }
};

} // namespace raptor

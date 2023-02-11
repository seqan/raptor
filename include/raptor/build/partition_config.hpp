// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides raptor::partition_config.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <bit>
#include <cstddef>
#include <cstdint>

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
        size_t const suffixes = next_power_of_four(partitions);
        size_t const suffixes_per_part = suffixes / partitions;
        mask = suffixes - 1;
        shift_value = std::countr_zero(suffixes_per_part);
    }

    size_t partitions{};
    size_t mask{};
    int shift_value{};

    // The number of suffixes is a power of four.
    // The number of partitions is a power of two.
    // The number of suffixes is greater than or equal to the number of partitions.
    // This means that the number of suffixes per partition is always a power of two.
    // Therefore, we can do a right shift instead of the division in:
    // (hash & mask) / suffixes_per_part == partition
    // The compiler cannot optimize the division to a right shift because suffixes_per_part is a runtime variable.
    constexpr size_t hash_partition(uint64_t const hash) const
    {
        return (hash & mask) >> shift_value;
    }

    static constexpr size_t next_power_of_four(size_t number)
    {
        if (number == 0ULL || number == 1ULL)
            return 1ULL; // GCOVR_EXCL_LINE

        --number;
        int const highest_set_bit_pos = std::bit_width(number);
        int const shift_amount = (highest_set_bit_pos + (highest_set_bit_pos & 1)) - 2;
        //                       (           Next multiple of two                )   4 has two zeros
        return 0b0100ULL << shift_amount;
    }
};

} // namespace raptor

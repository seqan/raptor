// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::formatted_bytes.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <string>

namespace raptor
{

[[nodiscard]] inline std::string formatted_bytes(size_t const bytes)
{
    assert(bytes > 0);

    size_t iterations{};
    size_t integer{bytes};

    while (integer >> 10u && iterations < 6u)
    {
        integer >>= 10u;
        ++iterations;
    }

    // While this is a bit more involved, we can avoid using floating point numbers.
    auto first_decimal_position = [&]()
    {
        assert(iterations > 0u);
        size_t decimal{bytes};
        decimal -= integer << (iterations * 10u);             // Substract bytes represented by integer, e.g. -5GiB
        decimal >>= (iterations - 1u) * 10u;                  // Shift to next smallest unit, e.g. 800MiB
        decimal = decimal * 1000u / 1024u;                    // Account for using decimal system, i.e. 800MiB != 0.8GiB
        size_t const diff{decimal - (decimal / 100u) * 100u}; // We want to round up to 1 decimal position
        bool const round_up{diff >= 50u};
        decimal += static_cast<size_t>(round_up * 100u) - diff;
        decimal /= 100u;
        return decimal;
    };

    auto formatted_string = [&]()
    {
        static constexpr int8_t int_to_char_offset{'0'}; // int 0 as char: char{0 + 48} = '0'
        size_t const decimal = iterations ? first_decimal_position() : 0u;
        assert(decimal <= 10u);

        if (!iterations) // No decimals for Bytes
            return std::to_string(integer);
        else if (decimal < 10u) // No need to round integer part
            return std::to_string(integer) + '.' + static_cast<char>(decimal + int_to_char_offset);
        else // Round integer part, e.g., 5.99 MiB should report 6.0 MiB
        {
            ++integer;
            // Check whether rounding results in a change of unit, e.g. 1023.99MiB to 1.0GiB
            if (integer >> 10u)
            {
                ++iterations;
                integer >>= 10u;
            }
            return std::to_string(integer) + ".0";
        }
    };

    std::string const formatted{formatted_string()};
    switch (iterations)
    {
    case 0:
        return "[Bytes]: " + formatted;
    case 1:
        return "[KiB]: " + formatted;
    case 2:
        return "[MiB]: " + formatted;
    case 3:
        return "[GiB]: " + formatted;
    case 4:
        return "[TiB]: " + formatted;
    case 5:
        return "[PiB]: " + formatted;
    default:
        return "[EiB]: " + formatted;
    }
}

} // namespace raptor

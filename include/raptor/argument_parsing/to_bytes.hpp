// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::to_bytes.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <algorithm>    // for min
#include <charconv>     // for from_chars
#include <cmath>        // for modf
#include <concepts>     // for same_as
#include <cstddef>      // for size_t
#include <memory>       // for allocator, addressof
#include <stdexcept>    // for invalid_argument, out_of_range
#include <string>       // for char_traits, operator+, basic_string, string
#include <string_view>  // for basic_string_view, operator==, string_view
#include <system_error> // for errc

#include <seqan3/std/charconv> // for from_chars

namespace raptor
{

namespace detail
{

template <typename number_t>
[[nodiscard]] size_t to_bytes_impl(std::string_view const input)
{
    number_t value{};
    auto [ptr, ec] = std::from_chars(input.data(), input.data() + input.size(), value);

    if (ec == std::errc::invalid_argument)
        throw std::invalid_argument{"Invalid number in \"" + std::string{input} + "\"."};
    else if (ec == std::errc::result_out_of_range)
        throw std::out_of_range{"Number in \"" + std::string{input} + "\" is out of range."};

    std::string_view const unit{ptr, input.end()};

    if (unit.empty())
    {
        if constexpr (std::same_as<number_t, double>)
        {
            double integral_part{};
            double fractional_part = std::modf(value, std::addressof(integral_part));
            if (fractional_part)
                throw std::invalid_argument{"Fractional bytes are not allowed."};
        }

        return static_cast<size_t>(value);
    }

    if (unit == "K")
        return static_cast<size_t>(value * 1000);
    else if (unit == "Ki")
        return static_cast<size_t>(value * 1024);
    else if (unit == "M")
        return static_cast<size_t>(value * 1000 * 1000);
    else if (unit == "Mi")
        return static_cast<size_t>(value * 1024 * 1024);
    else if (unit == "G")
        return static_cast<size_t>(value * 1000 * 1000 * 1000);
    else if (unit == "Gi")
        return static_cast<size_t>(value * 1024 * 1024 * 1024);
    else if (unit == "T")
        return static_cast<size_t>(value * 1000ULL * 1000ULL * 1000ULL * 1000ULL);
    else if (unit == "Ti")
        return static_cast<size_t>(value * 1024ULL * 1024ULL * 1024ULL * 1024ULL);
    else if (unit == "P")
        return static_cast<size_t>(value * 1000ULL * 1000ULL * 1000ULL * 1000ULL * 1000ULL);
    else if (unit == "Pi")
        return static_cast<size_t>(value * 1024ULL * 1024ULL * 1024ULL * 1024ULL * 1024ULL);
    else if (unit == "E")
        return static_cast<size_t>(value * 1000ULL * 1000ULL * 1000ULL * 1000ULL * 1000ULL * 1000ULL);
    else if (unit == "Ei")
        return static_cast<size_t>(value * 1024ULL * 1024ULL * 1024ULL * 1024ULL * 1024ULL * 1024ULL);

    throw std::invalid_argument{"Unknown unit in \"" + std::string{input} + "\"."};
}

} // namespace detail

[[nodiscard]] inline size_t to_bytes(std::string_view input)
{
    input.remove_prefix(std::min(input.find_first_not_of(' '), input.size()));
    if (auto const pos = input.find_last_not_of(' '); pos != std::string_view::npos)
        input.remove_suffix(input.size() - 1 - pos);

    if (input.empty())
        return 0;

    if (input.starts_with('-'))
        throw std::invalid_argument{"Negative numbers are not allowed."};

    if (input.contains('.'))
        return detail::to_bytes_impl<double>(input);
    else
        return detail::to_bytes_impl<size_t>(input);
}

} // namespace raptor

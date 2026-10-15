// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::json_writer.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <array>
#include <cassert>
#include <charconv>
#include <cmath>
#include <limits>
#include <ostream>

#include <raptor/inspect/json_writer.hpp>

namespace raptor
{

std::string json_quote(std::string_view const value)
{
    static constexpr std::string_view hex_digits{"0123456789abcdef"};

    std::string result{};
    result.reserve(value.size() + 2u);
    result += '"';

    for (char const character : value)
    {
        switch (character)
        {
        case '"':
            result += "\\\"";
            break;
        case '\\':
            result += "\\\\";
            break;
        case '\b':
            result += "\\b";
            break;
        case '\f':
            result += "\\f";
            break;
        case '\n':
            result += "\\n";
            break;
        case '\r':
            result += "\\r";
            break;
        case '\t':
            result += "\\t";
            break;
        default:
            // JSON forbids unescaped control characters. Those without a short escape are written as `\u00XX`.
            if (auto const byte = static_cast<unsigned char>(character); byte < 0x20u)
            {
                result += "\\u00";
                result += hex_digits[byte >> 4u];
                result += hex_digits[byte & 0x0Fu];
            }
            else
            {
                result += character;
            }
        }
    }

    result += '"';
    return result;
}

std::string shortest_double(double const value)
{
    // The longest shortest representation is 24 characters, e.g. -2.2250738585072014e-308.
    std::array<char, 32u> buffer{};
    auto const [end, error] = std::to_chars(buffer.data(), buffer.data() + buffer.size(), value);
    assert(error == std::errc{});
    return {buffer.data(), end};
}

json_writer & json_writer::key(std::string_view const name)
{
    // A key is only valid inside an object, and only once per value.
    assert(!open_containers.empty());
    assert(open_containers.back().is_object);
    assert(!has_key);

    write_separator(open_containers.back());
    stream << json_quote(name) << ": ";
    has_key = true;
    return *this;
}

json_writer & json_writer::begin_object(bool const inline_elements)
{
    open(true, inline_elements);
    return *this;
}

json_writer & json_writer::end_object()
{
    close(true);
    return *this;
}

json_writer & json_writer::begin_array(bool const inline_elements)
{
    open(false, inline_elements);
    return *this;
}

json_writer & json_writer::end_array()
{
    close(false);
    return *this;
}

json_writer & json_writer::string(std::string_view const value)
{
    prepare_value();
    stream << json_quote(value);
    return *this;
}

json_writer & json_writer::boolean(bool const value)
{
    prepare_value();
    stream << (value ? "true" : "false");
    return *this;
}

json_writer & json_writer::null()
{
    prepare_value();
    stream << "null";
    return *this;
}

json_writer & json_writer::write_unsigned(uint64_t const value)
{
    prepare_value();
    std::array<char, std::numeric_limits<uint64_t>::digits10 + 1> buffer{};
    auto const [end, error] = std::to_chars(buffer.data(), buffer.data() + buffer.size(), value);
    assert(error == std::errc{});
    stream.write(buffer.data(), end - buffer.data());
    return *this;
}

json_writer & json_writer::write_double(double const value)
{
    // JSON has no representation for NaN and infinity.
    if (!std::isfinite(value))
        return null();

    prepare_value();
    stream << shortest_double(value);
    return *this;
}

void json_writer::prepare_value()
{
    // The top-level value.
    if (open_containers.empty())
    {
        assert(!has_key);
        return;
    }

    // Inside an object, raptor::json_writer::key already wrote the separator.
    if (open_containers.back().is_object)
    {
        assert(has_key);
        has_key = false;
        return;
    }

    assert(!has_key);
    write_separator(open_containers.back());
}

void json_writer::write_separator(container & parent)
{
    if (parent.has_elements)
        stream << ',';

    if (parent.inline_elements)
    {
        if (parent.has_elements)
            stream << ' ';
    }
    else
    {
        stream << '\n';
        indent(open_containers.size());
    }

    parent.has_elements = true;
}

void json_writer::open(bool const is_object, bool const inline_elements)
{
    prepare_value();
    // Everything nested in an inline container is inline as well.
    bool const parent_is_inline = !open_containers.empty() && open_containers.back().inline_elements;
    open_containers.push_back({.is_object = is_object, .inline_elements = inline_elements || parent_is_inline});
    stream << (is_object ? '{' : '[');
}

void json_writer::close(bool const is_object)
{
    // Containers must be closed in the reverse order of opening, and a key must be followed by a value.
    assert(!open_containers.empty());
    assert(open_containers.back().is_object == is_object);
    assert(!has_key);

    container const closed = open_containers.back();
    open_containers.pop_back();

    if (closed.has_elements && !closed.inline_elements)
    {
        stream << '\n';
        indent(open_containers.size());
    }

    stream << (is_object ? '}' : ']');

    // Terminate the document.
    if (open_containers.empty())
        stream << '\n';
}

void json_writer::indent(size_t const level)
{
    for (size_t i = 0; i < level; ++i)
        stream << "  ";
}

} // namespace raptor

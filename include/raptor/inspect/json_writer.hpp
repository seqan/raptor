// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::json_writer.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <cstdint>
#include <iosfwd>
#include <string>
#include <string_view>
#include <vector>

namespace raptor
{

//!\brief Returns `value` as a JSON string literal, including the enclosing quotes.
std::string json_quote(std::string_view const value);

/*!\brief Returns the shortest representation of `value` that reads back as the same `double`.
 * \details Uses `std::to_chars`, e.g. `0.05` and `1`, not `0.050000` and `1.000000`. Not finite values yield `nan`,
 *          `inf` and `-inf`, which are not valid JSON; raptor::json_writer::number writes `null` instead.
 */
std::string shortest_double(double const value);

/*!\brief Writes JSON to a stream, indented by two spaces.
 * \details
 * The caller states the structure; the writer keeps track of nesting, separators and indentation. Inside an object,
 * every value must be preceded by raptor::json_writer::key; inside an array, it must not.
 *
 * A container opened with `inline_elements` is written on a single line, as are all containers nested inside it.
 * This keeps long lists of short entries, e.g. one per user bin, at one line per entry.
 *
 * Strings are written as given: bytes that are not valid UTF-8, e.g. in file names, are passed through.
 *
 * Raptor only writes JSON, a single small object for `raptor inspect`, which does not justify a library:
 * * simdjson is a compiled dependency, and its `string_builder` leaves the separators to the caller as well.
 * * The rapidjson bundled with cereal writes `NaN` and `Infinity`, which are not valid JSON, because cereal sets
 *   `CEREAL_RAPIDJSON_WRITE_DEFAULT_FLAGS` to `kWriteNanAndInfFlag`. Its `PrettyWriter` can only write arrays on a
 *   single line, not objects.
 *
 * \attention The output is only complete after every container has been closed. Correct nesting is checked by
 *            assertions only.
 */
class json_writer
{
public:
    json_writer() = delete;
    json_writer(json_writer const &) = delete;             // Would write to the same stream with a diverging state.
    json_writer & operator=(json_writer const &) = delete; // Would write to the same stream with a diverging state.
    json_writer(json_writer &&) = default;
    json_writer & operator=(json_writer &&) = delete; // reference member
    ~json_writer() = default;

    explicit json_writer(std::ostream & stream) : stream{stream}
    {}

    //!\brief Names the next value. Only valid inside an object.
    json_writer & key(std::string_view const name);

    json_writer & begin_object(bool const inline_elements = false);
    json_writer & end_object();
    json_writer & begin_array(bool const inline_elements = false);
    json_writer & end_array();

    json_writer & string(std::string_view const value);
    json_writer & boolean(bool const value);
    json_writer & null();

    template <std::unsigned_integral number_t>
    json_writer & number(number_t const value)
    {
        return write_unsigned(static_cast<uint64_t>(value));
    }

    //!\brief Writes `null` if `value` is not finite.
    template <std::floating_point number_t>
    json_writer & number(number_t const value)
    {
        return write_double(static_cast<double>(value));
    }

private:
    //!\brief The state of an open container.
    struct container
    {
        bool is_object{};
        bool inline_elements{};
        bool has_elements{};
    };

    std::ostream & stream;
    std::vector<container> open_containers{};
    //!\brief Whether raptor::json_writer::key has been called for the next value.
    bool has_key{false};

    //!\brief Writes what has to precede a value: a separator and indentation, or nothing after a key.
    void prepare_value();
    //!\brief Writes the comma and line break or space that precede an element of `parent`.
    void write_separator(container & parent);
    void open(bool const is_object, bool const inline_elements);
    void close(bool const is_object);
    void indent(size_t const level);
    json_writer & write_unsigned(uint64_t const value);
    json_writer & write_double(double const value);
};

} // namespace raptor

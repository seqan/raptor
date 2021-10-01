// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <seqan3/std/charconv>
#include <seqan3/std/ranges>

#include <raptor/build/hibf/parse_chopper_pack_line.hpp>
#include <raptor/string_view.hpp>

// LCOV_EXCL_START
#if defined(__GNUC__) && (__GNUC__ < 10)
namespace std::ranges
{

using ::ranges::common_view;
using ::ranges::split_view;

} // namespace std::ranges
#endif
// LCOV_EXCL_END

namespace raptor::hibf
{

chopper_pack_record parse_chopper_pack_line(std::string const & current_line)
{
    chopper_pack_record result{};

    // initialize parsing
    std::string_view const buffer{current_line};
    auto const buffer_end{buffer.end()};
    auto field_end = buffer.begin();
    while (field_end != buffer_end && *field_end != '\t') ++field_end;

    // parse filenames
    std::string_view const filenames{raptor::detail::string_view(buffer.begin(), field_end)};
    for (auto && filename : std::ranges::split_view{filenames, ';'})
    {
        auto const common_view = std::ranges::common_view{filename};
        result.filenames.emplace_back(common_view.begin(), common_view.end());
    }

    size_t tmp{};

    do // read bin_indices
    {
        ++field_end; // skip tab or ;
        field_end = std::from_chars(field_end, buffer_end, tmp).ptr;
        result.bin_indices.push_back(tmp);
    } while (field_end != buffer_end && *field_end != '\t');

    do // read number of technical bins
    {
        ++field_end; // skip tab or ;
        field_end = std::from_chars(field_end, buffer_end, tmp).ptr;
        result.number_of_bins.push_back(tmp);
    } while (field_end != buffer_end && *field_end != '\t');

    return result;
}

} // namespace raptor::hibf

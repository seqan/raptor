// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <seqan3/std/charconv>

#include <seqan3/utility/views/to.hpp>

#include <raptor/build/hibf/parse_chopper_pack_line.hpp>

namespace raptor::hibf
{

chopper_pack_record parse_chopper_pack_line(std::string const & current_line)
{
    chopper_pack_record result{};

    // initialize parsing
    char const * buffer = current_line.c_str();
    auto field_start = &buffer[0];
    auto field_end = &buffer[0];
    auto const buffer_end = field_start + current_line.size();
    while (field_end != buffer_end && *field_end != '\t') ++field_end;

    // parse filenames
    std::string filenames_str = std::string(field_start, field_end);
    for (auto && filename : filenames_str | std::views::split(';'))
        result.filenames.push_back((filename | seqan3::views::to<std::string>)); // replace TODO

    size_t tmp; // temporary size_t

    do // read bin_indices
    {
        ++field_end; // skip tab or ;
        auto res = std::from_chars(field_end, buffer_end, tmp);
        field_end = res.ptr;
        result.bin_indices.push_back(tmp);
    } while (field_end != buffer_end && *field_end != '\t');

    do // read number of technical bins
    {
        ++field_end; // skip tab or ;
        auto res = std::from_chars(field_end, buffer_end, tmp);
        field_end = res.ptr;
        result.number_of_bins.push_back(tmp);
    } while (field_end != buffer_end && *field_end != '\t');

    do // read estimated maximum technical bin size
    {
        ++field_end; // skip tab or ;
        auto res = std::from_chars(field_end, buffer_end, tmp);
        field_end = res.ptr;
        result.estimated_sizes.push_back(tmp);
    } while (field_end != buffer_end && *field_end != '\n');

    return result;
}

} // namespace raptor::hibf

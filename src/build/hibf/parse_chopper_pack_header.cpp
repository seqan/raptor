// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Implements raptor::hibf::parse_chopper_pack_header.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <lemon/list_graph.h> /// Must be first include.

#include <algorithm>
#include <cassert>
#include <seqan3/std/charconv>

#include <chopper/layout/layout.hpp>
#include <chopper/prefixes.hpp>

#include <raptor/build/hibf/parse_chopper_pack_header.hpp>

namespace raptor::hibf
{

std::pair<size_t, std::vector<chopper::layout::layout::max_bin>>
parse_chopper_pack_header(std::istream & chopper_pack_file)
{
    auto parse_bin_indices = [](std::string_view const & buffer)
    {
        std::vector<size_t> result;

        auto buffer_start = &buffer[0];
        auto const buffer_end = buffer_start + buffer.size();

        size_t tmp{};

        while (buffer_start < buffer_end)
        {
            buffer_start = std::from_chars(buffer_start, buffer_end, tmp).ptr;
            ++buffer_start; // skip ;
            result.push_back(tmp);
        }

        return result;
    }; // GCOVR_EXCL_LINE

    auto parse_first_bin = [](std::string_view const & buffer)
    {
        size_t tmp{};
        std::from_chars(&buffer[0], &buffer[0] + buffer.size(), tmp);
        return tmp;
    }; // GCOVR_EXCL_LINE

    std::string line;

    while (std::getline(chopper_pack_file, line) && line.size() >= 2
           && std::string_view{line}.substr(0, 1) == chopper::prefix::header
           && std::string_view{line}.substr(1, 1) == chopper::prefix::header_config)
        ; // skip config in header

    assert(line[0] == '#'); // we are reading header lines
    assert(line.substr(1, chopper::prefix::high_level.size()) == chopper::prefix::high_level);

    // parse High Level max bin index
    assert(line.substr(chopper::prefix::high_level.size() + 2, 11) == "max_bin_id:");
    std::string_view const hibf_max_bin_str{line.begin() + 27, line.end()};
    size_t const hibf_max_bin = parse_first_bin(hibf_max_bin_str);

    std::vector<chopper::layout::layout::max_bin> header_max_bins{};

    // first read and parse header records, in order to sort them before adding them to the graph
    while (std::getline(chopper_pack_file, line) && line.substr(0, 6) != "#FILES")
    {
        assert(line.substr(1, chopper::prefix::merged_bin.size()) == chopper::prefix::merged_bin);

        // parse header line
        std::string_view const indices_str{
            line.begin() + 1 /*#*/ + chopper::prefix::merged_bin.size() + 1 /*_*/,
            std::find(line.begin() + chopper::prefix::merged_bin.size() + 2, line.end(), ' ')};

        assert(line.substr(chopper::prefix::merged_bin.size() + indices_str.size() + 3, 11) == "max_bin_id:");
        std::string_view const max_id_str{line.begin() + chopper::prefix::merged_bin.size() + indices_str.size() + 14,
                                          line.end()};

        header_max_bins.emplace_back(parse_bin_indices(indices_str), parse_first_bin(max_id_str));
    }

    return {hibf_max_bin, header_max_bins};
}

} // namespace raptor::hibf

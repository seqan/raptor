// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <algorithm>
#include <cassert>
#include <seqan3/std/charconv>

#include <lemon/list_graph.h> /// Must be first include.

#include <chopper/prefixes.hpp>

#include <raptor/build/hibf/bin_prefixes.hpp>
#include <raptor/build/hibf/parse_chopper_pack_header.hpp>

namespace raptor::hibf
{

size_t parse_chopper_pack_header(lemon::ListDigraph & ibf_graph,
                                 lemon::ListDigraph::NodeMap<node_data> & node_map,
                                 std::istream & chopper_pack_file)
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

    assert(line[0] == '#');                                    // we are reading header lines
    assert(line.substr(1, hibf_prefix.size()) == hibf_prefix); // first line should always be High level IBF

    // parse High Level max bin index
    assert(line.substr(hibf_prefix.size() + 2, 11) == "max_bin_id:");
    std::string_view const hibf_max_bin_str{line.begin() + 27, line.end()};

    auto high_level_node = ibf_graph.addNode(); // high-level node = root node
    node_map.set(high_level_node, {0, parse_first_bin(hibf_max_bin_str), 0, lemon::INVALID, {}});

    std::vector<std::pair<std::vector<size_t>, size_t>> header_records{};

    // first read and parse header records, in order to sort them before adding them to the graph
    while (std::getline(chopper_pack_file, line) && line.substr(0, 6) != "#FILES")
    {
        assert(line.substr(1, merged_bin_prefix.size()) == merged_bin_prefix);

        // parse header line
        std::string_view const indices_str{line.begin() + 1 /*#*/ + merged_bin_prefix.size() + 1 /*_*/,
                                           std::find(line.begin() + merged_bin_prefix.size() + 2, line.end(), ' ')};

        assert(line.substr(merged_bin_prefix.size() + indices_str.size() + 3, 11) == "max_bin_id:");
        std::string_view const max_id_str{line.begin() + merged_bin_prefix.size() + indices_str.size() + 14,
                                          line.end()};

        header_records.emplace_back(parse_bin_indices(indices_str), parse_first_bin(max_id_str));
    }

    // sort records ascending by the number of bin indices (corresponds to the IBF levels)
    std::ranges::sort(header_records,
                      [](auto const & r, auto const & l)
                      {
                          return r.first.size() < l.first.size();
                      });

    for (auto const & [bin_indices, max_id] : header_records)
    {
        // we assume that the header lines are in the correct order
        // go down the tree until you find the matching parent
        auto it = bin_indices.begin();
        lemon::ListDigraph::Node current_node = high_level_node; // start at root

        while (it != (bin_indices.end() - 1))
        {
            for (lemon::ListDigraph::OutArcIt arc_it(ibf_graph, current_node); arc_it != lemon::INVALID; ++arc_it)
            {
                auto target = ibf_graph.target(arc_it);
                if (node_map[target].parent_bin_index == *it)
                {
                    current_node = target;
                    break;
                }
            }
            ++it;
        }

        auto new_node = ibf_graph.addNode();
        ibf_graph.addArc(current_node, new_node);
        node_map.set(new_node, {bin_indices.back(), max_id, 0, lemon::INVALID, {}});

        if (node_map[current_node].max_bin_index == bin_indices.back())
            node_map[current_node].favourite_child = new_node;
    }

    return header_records.size();
}

} // namespace raptor::hibf

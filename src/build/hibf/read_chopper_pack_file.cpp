// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Implements raptor::hibf::read_chopper_pack_file.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <lemon/list_graph.h> /// Must be first include.

#include <chopper/layout/layout.hpp>

#include <raptor/build/hibf/parse_chopper_pack_header.hpp>
#include <raptor/build/hibf/parse_chopper_pack_line.hpp>
#include <raptor/build/hibf/read_chopper_pack_file.hpp>
#include <raptor/build/hibf/update_content_node_data.hpp>
#include <raptor/build/hibf/update_header_node_data.hpp>

namespace raptor::hibf
{

void read_chopper_pack_file(build_data & data, std::string const & chopper_pack_filename)
{
    chopper::layout::layout hibf_layout{};

    std::ifstream chopper_pack_file{chopper_pack_filename};

    if (!chopper_pack_file.good() || !chopper_pack_file.is_open())
        throw std::logic_error{"Could not open file " + chopper_pack_filename + " for reading"}; // GCOVR_EXCL_LINE

    // parse header
    // -------------------------------------------------------------------------
    parse_chopper_pack_header(chopper_pack_file, hibf_layout);
    data.number_of_ibfs = hibf_layout.max_bins.size() + 1;
    // Add high level node
    auto high_level_node = data.ibf_graph.addNode(); // high-level node = root node
    data.node_map.set(high_level_node, {0, hibf_layout.top_level_max_bin_id, 0, lemon::INVALID, {}});

    update_header_node_data(std::move(hibf_layout.max_bins), data.ibf_graph, data.node_map);

    std::string current_line;
    while (std::getline(chopper_pack_file, current_line))
        hibf_layout.user_bins.emplace_back(parse_chopper_pack_line(current_line, data.filenames));

    data.number_of_user_bins = hibf_layout.user_bins.size();
    update_content_node_data(std::move(hibf_layout.user_bins), data.ibf_graph, data.node_map);

    data.resize();
}

} // namespace raptor::hibf

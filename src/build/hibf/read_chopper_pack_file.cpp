// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <lemon/list_graph.h> /// Must be first include.

#include <raptor/build/hibf/parse_chopper_pack_header.hpp>
#include <raptor/build/hibf/parse_chopper_pack_line.hpp>
#include <raptor/build/hibf/read_chopper_pack_file.hpp>

namespace raptor::hibf
{

template <seqan3::data_layout data_layout_mode>
void read_chopper_pack_file(build_data<data_layout_mode> & data, std::string const & chopper_pack_filename)
{
    std::ifstream chopper_pack_file{chopper_pack_filename};

    if (!chopper_pack_file.good() || !chopper_pack_file.is_open())
        throw std::logic_error{"Could not open file " + chopper_pack_filename + " for reading"}; // GCOVR_EXCL_LINE

    // parse header
    // -------------------------------------------------------------------------
    data.number_of_ibfs = parse_chopper_pack_header(data.ibf_graph, data.node_map, chopper_pack_file) + 1;

    // parse lines
    // -------------------------------------------------------------------------
    std::string current_line;
    size_t user_bins{};
    while (std::getline(chopper_pack_file, current_line))
    {
        ++user_bins;
        chopper_pack_record const && record = parse_chopper_pack_line(current_line);

        // go down the tree until you find the matching parent
        lemon::ListDigraph::Node current_node = data.ibf_graph.nodeFromId(0); // start at root

        for (size_t i = 0; i < record.bin_indices.size() - 1; ++i)
        {
            size_t const bin = record.bin_indices[i];
            size_t const num_tbs = record.number_of_bins[i];
            auto & current_data = data.node_map[current_node];

            // update number of technical bins in current_node-IBF
            current_data.number_of_technical_bins = std::max(current_data.number_of_technical_bins, bin + num_tbs);

#ifndef NDEBUG
            bool found_next_node{false}; // sanity check
#endif
            for (lemon::ListDigraph::OutArcIt arc_it(data.ibf_graph, current_node); arc_it != lemon::INVALID; ++arc_it)
            {
                auto target = data.ibf_graph.target(arc_it);
                if (data.node_map[target].parent_bin_index == bin)
                {
                    current_node = target;
#ifndef NDEBUG
                    found_next_node = true;
#endif
                    break;
                }
            }
            assert(found_next_node);
        }

        size_t const bin = record.bin_indices.back();
        size_t const num_tbs = record.number_of_bins.back();
        auto & current_data = data.node_map[current_node];

        // update number of technical bins in current_node-IBF
        current_data.number_of_technical_bins = std::max(current_data.number_of_technical_bins, bin + num_tbs);

        if (record.bin_indices.back() == current_data.max_bin_index)
            current_data.remaining_records.insert(current_data.remaining_records.begin(), record);
        else
            current_data.remaining_records.push_back(record);
    }

    data.number_of_user_bins = user_bins;
    data.resize();
}

template void read_chopper_pack_file<seqan3::data_layout::uncompressed>(build_data<seqan3::data_layout::uncompressed> &,
                                                                        std::string const &);

template void read_chopper_pack_file<seqan3::data_layout::compressed>(build_data<seqan3::data_layout::compressed> &,
                                                                      std::string const &);

} // namespace raptor::hibf

// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides raptor::hibf::parse_chopper_pack_header.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <iosfwd>

#include <raptor/build/hibf/node_data.hpp>

namespace raptor::hibf
{

size_t parse_chopper_pack_header(lemon::ListDigraph & ibf_graph,
                                 lemon::ListDigraph::NodeMap<node_data> & node_map,
                                 std::istream & chopper_pack_file);

} // namespace raptor::hibf

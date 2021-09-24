// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <string>

#include <raptor/build/hibf/build_data.hpp>

namespace raptor::hibf
{

// data needs to be passed from outside sind the graph in data cannot be moved
void read_chopper_pack_file(build_data & data, std::string const & chopper_pack_filename);

} // namespace raptor::hibf

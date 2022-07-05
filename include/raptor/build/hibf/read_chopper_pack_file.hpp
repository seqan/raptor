// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <string>

#include <raptor/build/hibf/build_data.hpp>

namespace raptor::hibf
{

template <seqan3::data_layout data_layout_mode>
void read_chopper_pack_file(build_data<data_layout_mode> & data, std::string const & chopper_pack_filename);

} // namespace raptor::hibf

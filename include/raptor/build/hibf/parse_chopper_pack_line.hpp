// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <string>

#include <raptor/build/hibf/chopper_pack_record.hpp>

namespace raptor::hibf
{

chopper_pack_record parse_chopper_pack_line(std::string const & current_line);

} // namespace raptor::hibf

// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <raptor/argument_parsing/build_arguments.hpp>
#include <raptor/argument_parsing/upgrade_arguments.hpp>

namespace raptor
{

void parse_bin_path(build_arguments & arguments);
void parse_bin_path(upgrade_arguments & arguments);

} // namespace raptor

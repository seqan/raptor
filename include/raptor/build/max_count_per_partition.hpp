// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <string>
#include <vector>

#include <raptor/argument_parsing/build_arguments.hpp>
#include <raptor/argument_parsing/upgrade_arguments.hpp>
#include <raptor/build/partition_config.hpp>

namespace raptor
{

std::vector<size_t> max_count_per_partition(partition_config const & cfg, build_arguments const & arguments);
std::vector<size_t> max_count_per_partition(partition_config const & cfg, upgrade_arguments const & arguments);

} // namespace raptor

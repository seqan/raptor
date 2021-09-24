// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <robin_hood.h>

#include <raptor/build/hibf/build_config.hpp>
#include <raptor/build/hibf/chopper_pack_record.hpp>

namespace raptor::hibf
{

robin_hood::unordered_flat_set<size_t> compute_kmers(build_config const & config, chopper_pack_record const & record);

void compute_kmers(robin_hood::unordered_flat_set<size_t> & kmers,
                   build_config const & config,
                   chopper_pack_record const & record);

} // namespace raptor::hibf

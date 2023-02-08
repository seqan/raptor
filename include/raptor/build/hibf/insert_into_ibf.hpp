// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <robin_hood.h>

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

#include <raptor/argument_parsing/build_arguments.hpp>
#include <raptor/build/hibf/chopper_pack_record.hpp>

namespace raptor::hibf
{

// automatically does naive splitting if number_of_bins > 1
void insert_into_ibf(robin_hood::unordered_flat_set<size_t> & parent_kmers,
                     robin_hood::unordered_flat_set<size_t> const & kmers,
                     size_t const number_of_bins,
                     size_t const bin_index,
                     seqan3::interleaved_bloom_filter<> & ibf,
                     bool is_root,
                     timer & fill_ibf_timer,
                     timer & merge_kmers_timer);

void insert_into_ibf(build_arguments const & arguments,
                     chopper_pack_record const & record,
                     seqan3::interleaved_bloom_filter<> & ibf);

} // namespace raptor::hibf

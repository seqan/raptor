// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides raptor::hibf::insert_into_ibf.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <robin_hood.h>

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

#include <chopper/layout/layout.hpp>

#include <raptor/build/hibf/build_data.hpp>

namespace raptor::hibf
{

// automatically does naive splitting if number_of_bins > 1
void insert_into_ibf(robin_hood::unordered_flat_set<uint64_t> const & kmers,
                     size_t const number_of_bins,
                     size_t const bin_index,
                     seqan3::interleaved_bloom_filter<> & ibf,
                     timer<concurrent::yes> & fill_ibf_timer);

void insert_into_ibf(build_data const & data,
                     chopper::layout::layout::user_bin const & record,
                     seqan3::interleaved_bloom_filter<> & ibf);

} // namespace raptor::hibf

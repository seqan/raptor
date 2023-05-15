// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides raptor::hibf::update_parent_kmers.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <robin_hood.h>

#include <raptor/argument_parsing/timer.hpp>

namespace raptor::hibf
{

void update_parent_kmers(robin_hood::unordered_flat_set<uint64_t> & parent_kmers,
                         robin_hood::unordered_flat_set<uint64_t> const & kmers,
                         timer<concurrent::yes> & merge_kmers_timer);

} // namespace raptor::hibf

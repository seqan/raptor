// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides raptor::hibf::input_base.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <robin_hood.h>

namespace raptor::hibf
{

// GCOVR_EXCL_START
class input_base
{
public:
    virtual void hash_into(size_t const, robin_hood::unordered_flat_set<size_t> &) const
    {
        throw std::runtime_error{"Access raptor::hibf::input_base::hash_into"};
    }
};
// GCOVR_EXCL_STOP

} // namespace raptor::hibf

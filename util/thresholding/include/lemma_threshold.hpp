// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/eseiler/minimizer_thresholds/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides lemma_threshold.
 */

#pragma once

#include <cinttypes>

struct lemma_threshold
{
    lemma_threshold() = default;
    lemma_threshold(lemma_threshold const &) = default;
    lemma_threshold & operator=(lemma_threshold const &) = default;
    lemma_threshold(lemma_threshold &&) = default;
    lemma_threshold & operator=(lemma_threshold &&) = default;
    ~lemma_threshold() = default;

    inline uint32_t threshold(uint64_t const text_length, uint64_t const window_size, uint16_t const errors) const
    {
        if ((errors + 1) * window_size > text_length + 1)
            return 0;

        return text_length + 1 - (errors + 1) * window_size;
    }
};

// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides strong types to protect interfaces.
 */

#pragma once

//!\brief Strong type for passing the window size.
struct window
{
    uint64_t v;
};
//!\brief Strong type for passing the kmer size.
struct kmer
{
    uint8_t v;
};
//!\brief Strong type for passing number of bins.
struct bins
{
    uint64_t v;
};
//!\brief Strong type for passing number of bits.
struct bits
{
    uint64_t v;
};
//!\brief Strong type for passing number of hash functions.
struct hashes
{
    uint64_t v;
};

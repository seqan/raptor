// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/eseiler/minimizer_thresholds/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides strong types to protect interfaces.
 */

#pragma once

//!\brief Strong type for passing the window size.
struct window { uint64_t v; };
//!\brief Strong type for passing the kmer size.
struct kmer { uint8_t v; };
//!\brief Strong type for passing number of bins.
struct bins { uint64_t v; };
//!\brief Strong type for passing number of bits.
struct bits { uint64_t v; };
//!\brief Strong type for passing number of hash functions.
struct hashes { uint64_t v; };

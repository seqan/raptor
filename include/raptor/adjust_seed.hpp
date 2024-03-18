// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::adjust_seed.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <cstdint>

namespace raptor
{

/*\brief Adjust the default seed such that it does not interfere with the IBF's hashing.
 *\param kmer_size The used k-mer size. For gapped shapes, this corresponds to the number of set bits (count()).
 *\details
 *
 * The hashing used with the IBF assumes that the input values are uniformly distributed.
 * However, we use a 64 bit seed, and unless the `kmer_size` is 32, not all 64 bits of the k-mers change.
 * Hence, we need to shift the seed to the right.
 *
 * For example, using 2-mers and a seed of length 8 bit, the values for the k-mers will only change for the last 4 bits:
 *
 * ```
 * seed = 1111'1011
 * kmer = 0000'XXXX
 * ```
 *
 * `seed XOR kmer` will then always have 4 leading ones.
 */
static inline constexpr uint64_t adjust_seed(uint8_t const kmer_size) noexcept
{
    return 0x8F3F73B5CF1C9ADEULL >> (64u - 2u * kmer_size);
}

} // namespace raptor

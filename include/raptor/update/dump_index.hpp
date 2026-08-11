// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::dump_index.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <raptor/index.hpp>

namespace raptor
{

/*!\brief Prints the structure of an HIBF to `std::cerr`. Debug helper.
 * \details
 * For each IBF, prints the user bin each of its technical bins holds, where `D` denotes a deleted bin and `M<id>`
 * a merged bin pointing at IBF `<id>`. Afterwards, prints the bin count and the occupancy of each IBF.
 */
void dump_index(raptor_index<index_structure::hibf> const & index);
//!\overload
void dump_index(seqan::hibf::hierarchical_interleaved_bloom_filter const & hibf);

} // namespace raptor

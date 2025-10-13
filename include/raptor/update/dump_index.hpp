// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::raptor_update.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <hibf/hierarchical_interleaved_bloom_filter.hpp> // for hierarchical_interleaved_bloom_filter

#include <raptor/index.hpp> // for hibf, raptor_index

namespace raptor
{

void dump_index(raptor_index<index_structure::hibf> const & index);
void dump_index(seqan::hibf::hierarchical_interleaved_bloom_filter const & hibf);

} // namespace raptor

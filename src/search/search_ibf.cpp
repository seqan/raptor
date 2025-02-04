// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::search_ibf.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <raptor/search/search_ibf.hpp>
#include <raptor/search/search_singular_ibf.hpp>

namespace raptor
{

void search_ibf(search_arguments const & arguments)
{
    auto index = raptor_index<index_structure::ibf>{};
    search_singular_ibf(arguments, std::move(index));
}

} // namespace raptor

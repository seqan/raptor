// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::search_hibf.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <raptor/search/search_hibf.hpp>
#include <raptor/search/search_singular_ibf.hpp>

namespace raptor
{

void search_hibf(search_arguments const & arguments)
{
    auto index = raptor_index<index_structure::hibf>{};
    search_singular_ibf(arguments, std::move(index));
}

} // namespace raptor

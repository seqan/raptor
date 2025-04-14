// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::search_hibf.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <ranges>  // for views
#include <utility> // for move

#include <raptor/argument_parsing/search_arguments.hpp> // for search_arguments
#include <raptor/index.hpp>                             // for raptor_index, hibf
#include <raptor/search/search_hibf.hpp>                // for search_hibf
#include <raptor/search/search_singular_ibf.hpp>        // for search_singular_ibf

namespace raptor
{

void search_hibf(search_arguments const & arguments)
{
    auto index = raptor_index<index_structure::hibf>{};
    search_singular_ibf(arguments, std::move(index));
}

} // namespace raptor

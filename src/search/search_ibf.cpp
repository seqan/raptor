// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::search_ibf.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <charconv> // for to_chars
#include <ranges>   // for views
#include <utility>  // for move

#include <raptor/argument_parsing/search_arguments.hpp> // for search_arguments
#include <raptor/index.hpp>                             // for raptor_index, ibf
#include <raptor/search/search_ibf.hpp>                 // for search_ibf
#include <raptor/search/search_singular_ibf.hpp>        // for search_singular_ibf

namespace raptor
{

void search_ibf(search_arguments const & arguments)
{
    auto index = raptor_index<index_structure::ibf>{};
    search_singular_ibf(arguments, std::move(index));
}

} // namespace raptor

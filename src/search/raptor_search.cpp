// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::raptor_search.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <hibf/misc/timer.hpp> // for concurrent_timer

#include <raptor/argument_parsing/search_arguments.hpp> // for search_arguments
#include <raptor/search/search_hibf.hpp>                // for search_hibf
#include <raptor/search/search_ibf.hpp>                 // for search_ibf
#include <raptor/search/search_partitioned_ibf.hpp>     // for search_partitioned_ibf

namespace raptor
{

void raptor_search(search_arguments const & arguments)
{
    arguments.complete_search_timer.start();

    if (arguments.is_hibf)
        search_hibf(arguments);
    else if (arguments.parts == 1u)
    {
#if RAPTOR_FPGA
        if (arguments.use_fpga)
            search_fpga(arguments);
        else
#endif
            search_ibf(arguments);
    }
    else
        search_partitioned_ibf(arguments);

    arguments.complete_search_timer.stop();

    return;
}

} // namespace raptor

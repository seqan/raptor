// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::raptor_build.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <raptor/build/build_hibf.hpp>
#include <raptor/build/build_ibf.hpp>
#include <raptor/build/raptor_build.hpp>

namespace raptor
{

void raptor_build(build_arguments const & arguments)
{
    if (arguments.is_hibf)
        build_hibf(arguments);
    else
        build_ibf(arguments);
}

} // namespace raptor

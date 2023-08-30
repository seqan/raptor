// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

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

// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <raptor/build/build_ibf.hpp>
#include <raptor/build/hibf/chopper_build.hpp>
#include <raptor/build/raptor_build.hpp>

namespace raptor
{

void raptor_build(build_arguments const & arguments)
{
    if (arguments.is_hibf)
        if (arguments.compressed)
            hibf::chopper_build<seqan3::data_layout::compressed>(arguments);
        else
            hibf::chopper_build<seqan3::data_layout::uncompressed>(arguments);
    else if (arguments.compressed)
        build_ibf<true>(arguments);
    else
        build_ibf<false>(arguments);
}

} // namespace raptor

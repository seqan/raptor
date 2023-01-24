// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <raptor/search/search_hibf.hpp>
#include <raptor/search/search_ibf.hpp>
#include <raptor/search/search_partitioned_ibf.hpp>

namespace raptor
{

void raptor_search(search_arguments const & arguments)
{
    if (arguments.is_hibf)
    {
        if (arguments.compressed)
            search_hibf<true>(arguments);
        else
            search_hibf<false>(arguments);
    }
    else if (arguments.parts == 1)
    {
        if (arguments.compressed)
            search_ibf<true>(arguments);
        else
            search_ibf<false>(arguments);
    }
    else
    {
        if (arguments.compressed)
            search_partitioned_ibf<true>(arguments);
        else
            search_partitioned_ibf<false>(arguments);
    }

    return;
}

} // namespace raptor

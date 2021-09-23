// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <raptor/search/compute_simple_model.hpp>
#include <raptor/search/minimiser_model.hpp>

namespace raptor
{

std::vector<size_t> compute_simple_model(search_arguments const & arguments)
{
    std::vector<size_t> precomp_thresholds;

    if (!arguments.threshold && !do_cerealisation_in(precomp_thresholds, arguments))
    {
        precomp_thresholds = precompute_threshold(arguments.pattern_size,
                                                  arguments.window_size,
                                                  arguments.shape,
                                                  arguments.errors,
                                                  arguments.tau);

        do_cerealisation_out(precomp_thresholds, arguments);
    }

    return precomp_thresholds;
}

} // namespace raptor

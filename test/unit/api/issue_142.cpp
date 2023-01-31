// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/search/kmer_index/shape.hpp>

#include <raptor/threshold/one_indirect_error_model.hpp>

TEST(issue, 142)
{
    std::vector<double> const destroyed_by_error =
        raptor::threshold::one_indirect_error_model(15u, 6u, seqan3::ungapped{4u});
}

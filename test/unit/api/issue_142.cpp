// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h> // for Test, TestInfo, TEST

#include <string> // for basic_string
#include <vector> // for vector

#include <seqan3/search/kmer_index/shape.hpp> // for ungapped

#include <raptor/threshold/one_indirect_error_model.hpp> // for one_indirect_error_model

TEST(issue, 142)
{
    std::vector<double> const destroyed_by_error =
        raptor::threshold::one_indirect_error_model(15u, 6u, seqan3::ungapped{4u});
}

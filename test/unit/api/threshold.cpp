// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <raptor/threshold/threshold.hpp>

static inline raptor::threshold::threshold_parameters const default_parameters{.window_size = 32,
                                                                               .shape = seqan3::ungapped{32u},
                                                                               .query_length = 250u,
                                                                               .p_max = 0.15,
                                                                               .fpr = 0.05,
                                                                               .tau = 0.9999};

TEST(kmer_lemma, no_error)
{
    raptor::threshold::threshold threshold{default_parameters};
    EXPECT_EQ(threshold.get(100u), 219u);
    EXPECT_EQ(threshold.get(100u), threshold.get(1000u));
}

TEST(kmer_lemma, with_error)
{
    auto threshold_params = default_parameters;
    threshold_params.errors = 1u;
    raptor::threshold::threshold threshold{threshold_params};
    EXPECT_EQ(threshold.get(100u), 187u);
    EXPECT_EQ(threshold.get(100u), threshold.get(1000u));
}

TEST(kmer_lemma, minimum)
{
    auto threshold_params = default_parameters;
    threshold_params.errors = 10u;
    raptor::threshold::threshold threshold{threshold_params};
    EXPECT_EQ(threshold.get(100u), 1u);
    EXPECT_EQ(threshold.get(100u), threshold.get(1000u));
}

TEST(minimiser, no_error)
{
    auto threshold_params = default_parameters;
    threshold_params.window_size = 50u;
    raptor::threshold::threshold threshold{threshold_params};
    size_t const min_count{11u};  // 219 / 19
    size_t const max_count{201u}; // 250 - 50 + 1

    for (size_t i = 0; i < min_count; ++i)
        EXPECT_EQ(threshold.get(i), threshold.get(min_count)) << i;

    EXPECT_EQ(threshold.get(max_count), threshold.get(max_count + 1u));
    EXPECT_EQ(threshold.get(max_count), threshold.get(1000u));

    for (size_t i = min_count; i < max_count; ++i)
    {
        if (i == 59u) // 63 vs 60
            continue;
        EXPECT_LE(threshold.get(i), threshold.get(i + 1u)) << i;
    }
}

TEST(minimiser, with_error)
{
    auto threshold_params = default_parameters;
    threshold_params.window_size = 50u;
    threshold_params.errors = 1u;
    raptor::threshold::threshold threshold{threshold_params};
    size_t const min_count{11u};  // 219 / 19
    size_t const max_count{201u}; // 250 - 50 + 1

    for (size_t i = 0; i < min_count; ++i)
        EXPECT_EQ(threshold.get(i), threshold.get(min_count)) << i;

    EXPECT_EQ(threshold.get(max_count), threshold.get(max_count + 1u));
    EXPECT_EQ(threshold.get(max_count), threshold.get(1000u));

    for (size_t i = min_count; i < max_count; ++i)
    {
        if (i == 59u) // 47 vs 44
            continue;
        EXPECT_LE(threshold.get(i), threshold.get(i + 1u)) << i;
    }
}

TEST(minimiser, minimum)
{
    auto threshold_params = default_parameters;
    threshold_params.window_size = 50u;
    threshold_params.query_length = 75u;
    threshold_params.errors = 4u;
    raptor::threshold::threshold threshold{threshold_params};
    size_t const min_count{2u};  // 44 / 19
    size_t const max_count{26u}; // 75 - 50 + 1

    for (size_t i = 0; i <= min_count; ++i)
        EXPECT_EQ(threshold.get(i), 1u) << i;

    EXPECT_EQ(threshold.get(max_count), threshold.get(max_count + 1u));
    EXPECT_EQ(threshold.get(max_count), threshold.get(1000u));

    for (size_t i = min_count; i < max_count; ++i)
        EXPECT_LE(threshold.get(i), threshold.get(i + 1u)) << i;
}

TEST(percentage, 100)
{
    auto threshold_params = default_parameters;
    threshold_params.percentage = 1.0;
    raptor::threshold::threshold threshold{threshold_params};
    EXPECT_EQ(threshold.get(100u), 100u);
    EXPECT_EQ(threshold.get(250u), 250u);
}

TEST(percentage, 50)
{
    auto threshold_params = default_parameters;
    threshold_params.percentage = 0.5;
    raptor::threshold::threshold threshold{threshold_params};
    EXPECT_EQ(threshold.get(100u), 50u);
    EXPECT_EQ(threshold.get(250u), 125u);
}

TEST(percentage, 0)
{
    auto threshold_params = default_parameters;
    threshold_params.percentage = 0.0;
    raptor::threshold::threshold threshold{threshold_params};
    EXPECT_EQ(threshold.get(100u), 1u);
    EXPECT_EQ(threshold.get(250u), 1u);
}

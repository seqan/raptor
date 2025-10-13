// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h> // for Test, Message, TestPartResult, CmpHelperEQ, CmpHelperEQ...

#include <cstddef> // for size_t
#include <string>  // for basic_string
#include <vector>  // for vector
#include <version> // for _LIBCPP_VERSION

#include <seqan3/search/kmer_index/shape.hpp> // for shape, ungapped

#include <raptor/threshold/threshold.hpp>            // for threshold
#include <raptor/threshold/threshold_parameters.hpp> // for threshold_parameters

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

TEST(minimiser, logspace_substract)
{
    raptor::threshold::threshold_parameters const parameters{.window_size = 15,
                                                             .shape = seqan3::ungapped{13u},
                                                             .query_length = 50u,
                                                             .errors = 2u,
                                                             .p_max = 0.15,
                                                             .fpr = 0.05,
                                                             .tau = 0.9999};

    raptor::threshold::threshold const threshold{parameters};
    // results for pseudorandom generators are implementation-defined
#ifdef _LIBCPP_VERSION
    std::vector<size_t> const expected{1u,  1u,  1u,  1u,  2u,  4u,  5u,  6u,  7u,  8u,  9u,  9u, 10u,
                                       11u, 12u, 13u, 14u, 15u, 16u, 17u, 18u, 19u, 20u, 22u, 23u};
#else
    std::vector<size_t> const expected{1u,  1u,  1u,  1u,  2u,  4u,  5u,  6u,  7u,  7u,  8u,  9u, 10u,
                                       11u, 12u, 13u, 14u, 15u, 16u, 17u, 18u, 19u, 20u, 22u, 23u};
#endif
    ASSERT_EQ(expected.size(), 37u - 12u);
    for (size_t i = 12u; i < 37u; ++i)
        EXPECT_EQ(threshold.get(i), expected[i - 12u]) << i;
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

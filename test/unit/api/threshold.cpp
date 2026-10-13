// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <cmath>
#include <utility>

#include <raptor/threshold/logspace.hpp>
#include <raptor/threshold/pascal_row.hpp>
#include <raptor/threshold/threshold.hpp>

static inline raptor::threshold::threshold_parameters const default_parameters{.window_size = 32,
                                                                               .shape = seqan3::ungapped{32u},
                                                                               .query_length = 250u,
                                                                               .p_max = 0.15,
                                                                               .fpr = 0.05,
                                                                               .tau = 0.9999};

TEST(pascal_row, exact)
{
    std::vector<double> const row = raptor::threshold::pascal_row(19u);
    std::vector<long long> const expected{1,     19,    171,   969,   3876,  11628, 27132, 50388, 75582, 92378,
                                          92378, 75582, 50388, 27132, 11628, 3876,  969,   171,   19,    1};

    ASSERT_EQ(row.size(), expected.size());
    for (size_t i = 0; i < row.size(); ++i)
        EXPECT_EQ(std::llround(std::exp(row[i])), expected[i]) << i;
}

TEST(pascal_row, lgamma)
{
    // The largest deviation for n <= 250 is about 6e-13.
    for (size_t n = 0; n <= 250u; ++n)
    {
        std::vector<double> const row = raptor::threshold::pascal_row(n);
        ASSERT_EQ(row.size(), n + 1u);

        double const log_n_factorial = std::lgamma(n + 1.0);
        for (size_t i = 0; i <= n; ++i)
        {
            double const expected = log_n_factorial - std::lgamma(i + 1.0) - std::lgamma(n - i + 1.0);
            EXPECT_NEAR(row[i], expected, 1e-9) << n << ' ' << i;
        }
    }
}

TEST(logspace, substract)
{
    // y / x above and below 1 / 2 selects the expm1 and the log1p branch, respectively.
    std::vector<std::pair<double, double>> const
        values{{1.0, 0.2}, {0.8, 0.2}, {3.0, 1.0}, {1e-3, 1e-9}, {0.5, 0.49}, {2.0, 1.5}, {1.0, 0.6}};
    for (auto const & [x, y] : values)
        EXPECT_NEAR(raptor::logspace::substract(std::log(x), std::log(y)), std::log(x - y), 1e-12) << x << ' ' << y;

    EXPECT_DOUBLE_EQ(raptor::logspace::substract(std::log(0.8), raptor::logspace::negative_inf), std::log(0.8));
    EXPECT_EQ(raptor::logspace::substract(std::log(0.8), std::log(0.8)), raptor::logspace::negative_inf);
}

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
        if (i == 59u) // 44 vs 41
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
    // Results for pseudorandom generators are implementation-defined, but libc++ and libstdc++ agree here.
    std::vector<size_t> const expected{1u, 1u, 1u, 1u, 1u, 2u, 2u, 2u, 2u, 2u, 2u, 2u, 2u,
                                       2u, 2u, 2u, 3u, 4u, 4u, 5u, 6u, 7u, 8u, 9u, 10u};
    ASSERT_EQ(expected.size(), 37u - 12u);
    for (size_t i = 12u; i < 37u; ++i)
        EXPECT_EQ(threshold.get(i), expected[i - 12u]) << i;
}

TEST(minimiser, precompute_threshold)
{
    raptor::threshold::threshold_parameters const parameters{.window_size = 23,
                                                             .shape = seqan3::ungapped{19u},
                                                             .query_length = 50u,
                                                             .errors = 1u,
                                                             .p_max = 0.15,
                                                             .fpr = 0.05,
                                                             .tau = 0.9999};

    // One entry per number of minimisers, from 6 (32 / 5) to 28 (50 - 23 + 1).
    // Results for pseudorandom generators are implementation-defined, but libc++ and libstdc++ agree here.
    std::vector<size_t> const expected{0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u,
                                       0u, 0u, 1u, 2u, 2u, 3u, 4u, 4u, 5u, 6u, 7u};
    EXPECT_EQ(raptor::threshold::precompute_threshold(parameters), expected);
}

TEST(minimiser, precompute_correction)
{
    raptor::threshold::threshold_parameters const parameters{.window_size = 40,
                                                             .shape = seqan3::ungapped{20u},
                                                             .query_length = 99u,
                                                             .p_max = 0.15,
                                                             .fpr = 0.05,
                                                             .tau = 0.9999};

    // One entry per number of minimisers, from 3 (80 / 21) to 60 (99 - 40 + 1).
    std::vector<size_t> const expected{0u, 1u, 1u, 1u, 1u, 1u, 1u, 1u, 1u, 1u, 1u, 1u, 1u, 1u, 2u, 2u, 2u, 2u, 2u, 2u,
                                       2u, 2u, 2u, 2u, 2u, 2u, 2u, 2u, 2u, 2u, 2u, 3u, 3u, 3u, 3u, 3u, 3u, 3u, 3u, 3u,
                                       3u, 3u, 3u, 3u, 3u, 3u, 3u, 3u, 3u, 3u, 3u, 4u, 4u, 4u, 4u, 4u, 4u, 0u};
    EXPECT_EQ(raptor::threshold::precompute_correction(parameters), expected);
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

// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::threshold::threshold_parameters.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <cstdint>    // for uint32_t, uint64_t, uint8_t
#include <filesystem> // for path
#include <limits>     // for numeric_limits

#include <seqan3/search/kmer_index/shape.hpp> // for shape

namespace raptor::threshold
{

struct threshold_parameters
{
    // Basic.
    uint32_t window_size{};
    seqan3::shape shape{};
    uint64_t query_length{};

    // Threshold.
    uint8_t errors{};                                            // threshold_kinds::(probabilistic|lemma)
    double percentage{std::numeric_limits<double>::quiet_NaN()}; // threshold_kinds::percentage
    double p_max{};                                              // threshold_kinds::probabilistic
    double fpr{};                                                // threshold_kinds::probabilistic
    double tau{};                                                // threshold_kinds::probabilistic

    // Cache results.
    bool cache_thresholds{};
    std::filesystem::path output_directory{};
};

} // namespace raptor::threshold

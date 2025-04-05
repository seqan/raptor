// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::threshold::threshold.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <raptor/threshold/precompute_correction.hpp>
#include <raptor/threshold/precompute_threshold.hpp>

namespace raptor::threshold
{

class threshold
{
public:
    threshold() = default;
    threshold(threshold const &) = default;
    threshold & operator=(threshold const &) = default;
    threshold(threshold &&) = default;
    threshold & operator=(threshold &&) = default;
    ~threshold() = default;

    threshold(threshold_parameters const & arguments);

    size_t get(size_t const minimiser_count) const noexcept;

protected:
    enum class threshold_kinds : uint8_t
    {
        probabilistic,
        lemma,
        percentage
    };

    threshold_kinds threshold_kind{threshold_kinds::probabilistic};
    std::vector<size_t> precomp_correction{};
    std::vector<size_t> precomp_thresholds{};
    size_t kmer_lemma{};
    size_t minimal_number_of_minimizers{};
    size_t maximal_number_of_minimizers{};
    double threshold_percentage{};
};

} // namespace raptor::threshold

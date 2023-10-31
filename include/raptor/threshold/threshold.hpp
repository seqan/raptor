// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

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

private:
    enum class threshold_kinds
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

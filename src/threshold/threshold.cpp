// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::threshold::threshold.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <hibf/platform.hpp>

#include <raptor/threshold/threshold.hpp>

namespace raptor::threshold
{

threshold::threshold(threshold_parameters const & arguments)
{
    size_t const kmer_size{arguments.shape.size()};
    size_t const kmers_per_window = arguments.window_size - kmer_size + 1;

    if (!std::isnan(arguments.percentage))
    {
        threshold_kind = threshold_kinds::percentage;
        threshold_percentage = arguments.percentage;
    }
    else if (kmers_per_window == 1u)
    {
        threshold_kind = threshold_kinds::lemma;
        size_t const kmer_lemma_minuend = arguments.query_length + 1u;
        size_t const kmer_lemma_subtrahend = (arguments.errors + 1u) * kmer_size;
        kmer_lemma = kmer_lemma_minuend > kmer_lemma_subtrahend ? kmer_lemma_minuend - kmer_lemma_subtrahend : 1u;
    }
    else
    {
        threshold_kind = threshold_kinds::probabilistic;
        size_t const kmers_per_pattern = arguments.query_length - kmer_size + 1;
        minimal_number_of_minimizers = kmers_per_pattern / kmers_per_window;
        maximal_number_of_minimizers = arguments.query_length - arguments.window_size + 1;
        precomp_correction = precompute_correction(arguments);
        precomp_thresholds = precompute_threshold(arguments);
    }
}

size_t threshold::get(size_t const minimiser_count) const noexcept
{
    switch (threshold_kind)
    {
    case threshold_kinds::lemma:
        return kmer_lemma;
    case threshold_kinds::percentage:
        return std::max<size_t>(1u, minimiser_count * threshold_percentage);
    case threshold_kinds::probabilistic:
    {
        size_t const index = std::clamp(minimiser_count, minimal_number_of_minimizers, maximal_number_of_minimizers)
                           - minimal_number_of_minimizers;
        return std::max<size_t>(1u, precomp_thresholds[index] + precomp_correction[index]);
    }
    default: // GCOVR_EXCL_LINE
        HIBF_UNREACHABLE;
    }
}

} // namespace raptor::threshold

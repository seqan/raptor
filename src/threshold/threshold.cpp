// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::threshold::threshold.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <algorithm>

#include <hibf/misc/unreachable.hpp>

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
        // Without errors, every minimiser is required. There is nothing to precompute.
        if (errors = arguments.errors; errors == 0u)
            return;
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
        // Also for counts outside [minimal_number_of_minimizers, maximal_number_of_minimizers].
        if (errors == 0u)
            return std::max<size_t>(1u, minimiser_count);

        size_t const number_of_minimisers =
            std::clamp(minimiser_count, minimal_number_of_minimizers, maximal_number_of_minimizers);
        size_t const index = number_of_minimisers - minimal_number_of_minimizers;
        // The correction accounts for false positives among the minimisers affected by errors. It must not require more
        // minimisers than there are.
        size_t const corrected_threshold{precomp_thresholds[index] + precomp_correction[index]};
        return std::max<size_t>(1u, std::min(number_of_minimisers, corrected_threshold));
    }
    default: // GCOVR_EXCL_LINE
        seqan::hibf::unreachable();
    }
}

} // namespace raptor::threshold

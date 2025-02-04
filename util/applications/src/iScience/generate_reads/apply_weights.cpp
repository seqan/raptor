// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <cmath>
#include <concepts>
#include <cstddef>
#include <numeric>
#include <vector>

namespace raptor::util::generate_reads
{

// Input
// vec = arguments.number_of_reads_per_bin
// vec was initialised with the weights of the bins
// Output
// vec is modified in place to contain the weights normalized to the number of reads
void apply_weights(std::vector<size_t> & vec, size_t const number_of_reads)
{
    // std::reduce is the same as std::accumulate, except out of order execution is allowed
    // std::reduce also infers the type that should be used for the initial value (here size_t)
    auto const sum_of_weights = std::reduce(vec.begin(), vec.end());
    static_assert(std::same_as<decltype(sum_of_weights), size_t const>); // Distrust

    double const scaling_factor = static_cast<double>(number_of_reads) / sum_of_weights;

    for (size_t & weight : vec)
        weight = std::round(weight * scaling_factor);
}

} // namespace raptor::util::generate_reads

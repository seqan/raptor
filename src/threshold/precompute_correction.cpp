// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::threshold::precompute_correction.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <cassert>    // for assert
#include <cmath>      // for log, isnan
#include <cstddef>    // for size_t
#include <cstdint>    // for uint8_t
#include <filesystem> // for path, operator/, exists
#include <fstream>    // for char_traits, basic_ostream, operator<<, basic_ifstream
#include <sstream>    // for basic_stringstream
#include <string>     // for basic_string, string, allocator
#include <vector>     // for vector

#include <cereal/archives/binary.hpp> // for BinaryInputArchive, BinaryOutputArchive
#include <cereal/specialize.hpp>      // for specialization

#include <seqan3/search/kmer_index/shape.hpp> // for shape

#include <raptor/threshold/pascal_row.hpp>            // for pascal_row
#include <raptor/threshold/precompute_correction.hpp> // for precompute_correction
#include <raptor/threshold/threshold_parameters.hpp>  // for threshold_parameters

namespace raptor::threshold
{

[[nodiscard]] std::string const correction_filename(threshold_parameters const & arguments)
{
    std::stringstream stream{};
    stream << "correction_" << std::hex << arguments.query_length << '_' << arguments.window_size << '_'
           << arguments.shape.to_ulong() << '_' << arguments.p_max << '_' << arguments.fpr << ".bin";
    std::string result = stream.str();
    if (auto it = result.find("0."); it != std::string::npos)
        result.replace(it, 2, "");
    if (auto it = result.find("0."); it != std::string::npos)
        result.replace(it, 2, "");
    return result;
}

void write_correction(std::vector<size_t> const & vec, threshold_parameters const & arguments)
{
    if (!arguments.cache_thresholds)
        return;

    std::filesystem::path filename = arguments.output_directory / correction_filename(arguments);
    std::ofstream os{filename, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(vec);
}

bool read_correction(std::vector<size_t> & vec, threshold_parameters const & arguments)
{
    std::filesystem::path filename = arguments.output_directory / correction_filename(arguments);
    if (!arguments.cache_thresholds || !std::filesystem::exists(filename))
        return false;

    std::ifstream is{filename, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    iarchive(vec);
    return true;
}

[[nodiscard]] std::vector<size_t> precompute_correction(threshold_parameters const & arguments)
{
    uint8_t const kmer_size{arguments.shape.size()};
    assert(arguments.window_size != kmer_size); // Use k-mer lemma.
    assert(std::isnan(arguments.percentage));   // Use percentage.

    std::vector<size_t> correction;

    if (read_correction(correction, arguments))
        return correction;

    double const fpr{std::log(arguments.fpr)};
    double const inv_fpr{std::log(1.0 - arguments.fpr)};
    double const log_p_max{std::log(arguments.p_max)};
    size_t const kmers_per_window{arguments.window_size - kmer_size + 1};
    size_t const kmers_per_pattern{arguments.query_length - kmer_size + 1};
    size_t const minimal_number_of_minimisers{kmers_per_pattern / kmers_per_window};
    size_t const maximal_number_of_minimisers{arguments.query_length - arguments.window_size + 1};

    correction.reserve(maximal_number_of_minimisers - minimal_number_of_minimisers + 1);

    auto binom = [&fpr, &inv_fpr](std::vector<double> const & binom_coeff,
                                  size_t const number_of_minimisers,
                                  size_t const number_of_fp)
    {
        return binom_coeff[number_of_fp] + number_of_fp * fpr + (number_of_minimisers - number_of_fp) * inv_fpr;
    };

    // Iterate over the possible number of minimisers.
    for (size_t number_of_minimisers = minimal_number_of_minimisers;
         number_of_minimisers <= maximal_number_of_minimisers;
         ++number_of_minimisers)
    {
        size_t number_of_fp{1u};
        std::vector<double> const binom_coeff{pascal_row(number_of_minimisers)};
        // How many FPs to expect for a given fpr and number of minimisers?
        // The probability of seeing this many FP must be below p_max.
        while (binom(binom_coeff, number_of_minimisers, number_of_fp) >= log_p_max)
            ++number_of_fp; // GCOVR_EXCL_LINE

        correction.push_back(number_of_fp - 1);
    }
    assert(correction.size() != 0);

    write_correction(correction, arguments);

    return correction;
}

} // namespace raptor::threshold

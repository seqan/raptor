// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Implements raptor::threshold::precompute_threshold.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <fstream>

#include <cereal/types/vector.hpp>

#include <raptor/threshold/logspace.hpp>
#include <raptor/threshold/multiple_error_model.hpp>
#include <raptor/threshold/one_error_model.hpp>
#include <raptor/threshold/one_indirect_error_model.hpp>
#include <raptor/threshold/precompute_threshold.hpp>

namespace raptor::threshold
{

[[nodiscard]] std::string const threshold_filename(threshold_parameters const & arguments)
{
    std::stringstream stream{};
    stream << "threshold_" << std::hex << arguments.query_length << '_' << arguments.window_size << '_'
           << arguments.shape.to_ulong() << '_' << static_cast<uint16_t>(arguments.errors) << '_' << arguments.tau
           << ".bin";
    std::string result = stream.str();
    if (auto it = result.find("0."); it != std::string::npos)
        result.replace(it, 2, "");
    return result;
}

void write_thresholds(std::vector<size_t> const & vec, threshold_parameters const & arguments)
{
    if (!arguments.cache_thresholds)
        return;

    std::filesystem::path filename = arguments.output_directory / threshold_filename(arguments);
    std::ofstream os{filename, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(vec);
}

bool read_thresholds(std::vector<size_t> & vec, threshold_parameters const & arguments)
{
    std::filesystem::path filename = arguments.output_directory / threshold_filename(arguments);
    if (!arguments.cache_thresholds || !std::filesystem::exists(filename))
        return false;

    std::ifstream is{filename, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    iarchive(vec);
    return true;
}

[[nodiscard]] std::vector<size_t> precompute_threshold(threshold_parameters const & arguments)
{
    uint8_t const kmer_size{arguments.shape.size()};
    assert(arguments.window_size != kmer_size); // Use k-mer lemma.
    assert(std::isnan(arguments.percentage));   // Use percentage.

    std::vector<size_t> thresholds;

    if (read_thresholds(thresholds, arguments))
        return thresholds;

    double const log_tau{std::log(arguments.tau)};
    size_t const kmers_per_window{arguments.window_size - kmer_size + 1};
    size_t const kmers_per_pattern{arguments.query_length - kmer_size + 1};
    size_t const minimal_number_of_minimisers{kmers_per_pattern / kmers_per_window};
    size_t const maximal_number_of_minimisers{arguments.query_length - arguments.window_size + 1};

    thresholds.reserve(maximal_number_of_minimisers - minimal_number_of_minimisers + 1);

    // Probability that i minimisers are indirectly affected by one error.
    std::vector<double> const affected_by_one_error_indirectly_prob{
        one_indirect_error_model(arguments.query_length, arguments.window_size, arguments.shape)};

    // Iterate over the possible number of minimisers.
    for (size_t number_of_minimisers = minimal_number_of_minimisers;
         number_of_minimisers <= maximal_number_of_minimisers;
         ++number_of_minimisers)
    {
        // Probability that a minimiser starts at index i. Uniform => Number of minimisers / possible indices.
        double const uniform_start_index_prob{std::log(number_of_minimisers) - std::log(kmers_per_pattern)};

        // Probability that i minimisers are affected by one error (directly or indirectly).
        std::vector<double> const affected_by_one_error_prob{
            one_error_model(kmer_size, uniform_start_index_prob, affected_by_one_error_indirectly_prob)};

        // Probability that i minimisers are affected e errors.
        std::vector<double> const affected_by_e_errors_prob{
            multiple_error_model(number_of_minimisers, arguments.errors, affected_by_one_error_prob)};

        // Max number of affected minimisers as predicted by `multiple_error_model`. Used for a check when adding
        // the probabilities.
        // This check is not strictly necessary, but in case of floating point number inaccuracies, it prevents
        // adding all probabilities in `affected_by_e_errors_prob`.
        // While `affected_by_e_errors_prob` computes all probabilities according to a theoretical worst case,
        // in practice, there are probabilites of 0 starting at a certain number of affected minimisers.
        size_t const max_affected =
            std::ranges::find(affected_by_e_errors_prob, logspace::negative_inf) - affected_by_e_errors_prob.begin();

        // The fraction of covered cases.
        double cumulative_prob{affected_by_e_errors_prob[0]};
        // How many minimisers are affected at most...
        size_t affected_minimisers{};
        // such that threshold holds with a probability of at least (1 - tau)?
        while (cumulative_prob < log_tau && affected_minimisers < max_affected)
            cumulative_prob = logspace::add(cumulative_prob, affected_by_e_errors_prob[++affected_minimisers]);

        assert(affected_minimisers <= number_of_minimisers);
        // Hence, there are at least this many left unaffected (threshold).
        thresholds.push_back(number_of_minimisers - affected_minimisers);
    }
    assert(thresholds.size() == maximal_number_of_minimisers - minimal_number_of_minimisers + 1);

    write_thresholds(thresholds, arguments);

    return thresholds;
}

} // namespace raptor::threshold

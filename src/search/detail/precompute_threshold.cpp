// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides stuff.
 */

#include <fstream>

#include <cereal/types/vector.hpp>

#include <raptor/search/detail/logspace.hpp>
#include <raptor/search/detail/multiple_error_model.hpp>
#include <raptor/search/detail/one_error_model.hpp>
#include <raptor/search/detail/one_indirect_error_model.hpp>
#include <raptor/search/detail/precompute_threshold.hpp>

namespace raptor::detail
{

[[nodiscard]] std::string const threshold_filename(search_arguments const & arguments)
{
    std::stringstream stream{};
    stream << "threshold_"
           << std::hex
           << arguments.pattern_size
           << '_'
           << arguments.window_size
           << '_'
           << arguments.shape.to_ulong()
           << '_'
           << static_cast<uint16_t>(arguments.errors)
           << '_'
           << arguments.tau
           << ".bin";
    std::string result = stream.str();
    if (auto it = result.find("0."); it != std::string::npos)
        result.replace(it, 2, "");
    return result;
}

void write_thresholds(std::vector<size_t> const & vec, search_arguments const & arguments)
{
    if (!arguments.cache_thresholds)
        return;

    std::filesystem::path filename = arguments.index_file.parent_path() / threshold_filename(arguments);
    std::ofstream os{filename, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(vec);
}

bool read_thresholds(std::vector<size_t> & vec, search_arguments const & arguments)
{
    std::filesystem::path filename = arguments.index_file.parent_path() / threshold_filename(arguments);
    if (!arguments.cache_thresholds || !std::filesystem::exists(filename))
        return false;

    std::ifstream is{filename, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    iarchive(vec);
    return true;
}

[[nodiscard]] std::vector<size_t> precompute_threshold(search_arguments const & arguments)
{
    uint8_t const kmer_size{arguments.shape.size()};
    assert(arguments.window_size != kmer_size); // Use k-mer lemma.
    assert(!arguments.treshold_was_set); // Use percentage.

    std::vector<size_t> thresholds;

    if (read_thresholds(thresholds, arguments))
        return thresholds;

    double const log_tau{std::log(arguments.tau)};
    size_t const kmers_per_window{arguments.window_size - kmer_size + 1};
    size_t const kmers_per_pattern{arguments.pattern_size - kmer_size + 1};
    size_t const minimal_number_of_minimisers{kmers_per_pattern / kmers_per_window};
    size_t const maximal_number_of_minimisers{arguments.pattern_size - arguments.window_size + 1};

    thresholds.reserve(maximal_number_of_minimisers - minimal_number_of_minimisers + 1);

    // Probability that i minimisers are indirectly affected by one error.
    std::vector<double> const affected_by_one_error_indirectly_prob{
        detail::one_indirect_error_model(arguments.pattern_size,
                                         arguments.window_size,
                                         arguments.shape)
    };

    // Iterate over the possible number of minimisers.
    for (size_t number_of_minimisers = minimal_number_of_minimisers;
         number_of_minimisers <= maximal_number_of_minimisers;
         ++number_of_minimisers)
    {
        // Probability that a minimiser starts at index i. Uniform => Number of minimisers / possible indices.
        double const uniform_start_index_prob{std::log(number_of_minimisers) - std::log(kmers_per_pattern)};

        // Probability that i minimisers are affected by one error (directly or indirectly).
        std::vector<double> const affected_by_one_error_prob{
            detail::one_error_model(kmer_size,
                                    uniform_start_index_prob,
                                    affected_by_one_error_indirectly_prob)
        };

        // Probability that i minimisers are affected e errors.
        std::vector<double> const affected_by_e_errors_prob{
            detail::multiple_error_model(number_of_minimisers,
                                         arguments.errors,
                                         affected_by_one_error_prob)
        };

        // The fraction of covered cases.
        double cumulative_prob{affected_by_e_errors_prob[0]};
        // How many minimisers are affected at most...
        size_t affected_minimisers{};
        // such that threshold holds with a probability of at least (1 - tau)?
        while (cumulative_prob < log_tau)
            cumulative_prob = logspace::add(cumulative_prob, affected_by_e_errors_prob[++affected_minimisers]);

        assert(affected_minimisers <= number_of_minimisers);
        // Hence, there are at least this many left unaffected (threshold).
        thresholds.push_back(number_of_minimisers - affected_minimisers);
    }
    assert(thresholds.size() == maximal_number_of_minimisers - minimal_number_of_minimisers + 1);

    write_thresholds(thresholds, arguments);

    return thresholds;
}

} // namespace raptor::detail

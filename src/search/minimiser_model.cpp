// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/eseiler/minimizer_thresholds/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides stuff.
 */

#include <numeric>

#include <cereal/types/vector.hpp>

#include <raptor/search/minimiser_model.hpp>

namespace raptor
{

std::vector<size_t> precompute_threshold(size_t const pattern_size,
                                         size_t const window_size,
                                         uint8_t const kmer_size,
                                         size_t const errors,
                                         double const tau)
{
    if (window_size == kmer_size)
        return {pattern_size + 1 > (errors + 1) * kmer_size ? pattern_size + 1 - (errors + 1) * kmer_size : 0};

    std::vector<size_t> thresholds;
    size_t const kmers_per_window = window_size - kmer_size + 1;
    size_t const kmers_per_pattern = pattern_size - kmer_size + 1;

    size_t const minimal_number_of_minimizers = std::ceil(kmers_per_pattern / static_cast<double>(kmers_per_window));
    size_t const maximal_number_of_minimizers = pattern_size - window_size + 1;

    std::vector<double> indirect_errors;
    indirect_errors = detail::destroyed_indirectly_by_error(pattern_size, window_size, kmer_size);

    // Iterate over the possible number of minimizers
    for (size_t number_of_minimizers = minimal_number_of_minimizers; number_of_minimizers <= maximal_number_of_minimizers; ++number_of_minimizers)
    {
        std::vector<double> proba_x(kmers_per_pattern, number_of_minimizers / static_cast<double>(kmers_per_pattern));

        auto [p_mean, proba] = detail::simple_model(kmer_size, proba_x, indirect_errors);
        (void) p_mean;

        std::vector<double> proba_error(number_of_minimizers, 0);
        for (size_t i = 0; i < number_of_minimizers; ++i)
            proba_error[i] = detail::enumerate_all_errors(i, errors, proba);

        double sum = std::accumulate(proba_error.begin(), proba_error.end(), 0.0);
        for (auto & x : proba_error)
            x /= sum;

        double n =0;
        for (size_t i = 0; i < number_of_minimizers; ++i)
        {
            n += proba_error[i];

            if (n >= tau)
            {
                thresholds.push_back(number_of_minimizers - i);
                break;
            }
        }
    }
    assert(thresholds.size() != 0);
    return thresholds;
}

void do_cerealisation_out(std::vector<size_t> const & vec, search_arguments const & arguments)
{
    std::filesystem::path filename = arguments.ibf_file.parent_path() / ("binary_p" + std::to_string(arguments.pattern_size) +
                                                                         "_w" + std::to_string(arguments.window_size) +
                                                                         "_k" + std::to_string(arguments.kmer_size) +
                                                                         "_e" + std::to_string(arguments.errors) +
                                                                         "_tau" + std::to_string(arguments.tau));
    std::ofstream os{filename, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(vec);
}

bool do_cerealisation_in(std::vector<size_t> & vec, search_arguments const & arguments)
{
    std::filesystem::path filename = arguments.ibf_file.parent_path() / ("binary_p" + std::to_string(arguments.pattern_size) +
                                                                         "_w" + std::to_string(arguments.window_size) +
                                                                         "_k" + std::to_string(arguments.kmer_size) +
                                                                         "_e" + std::to_string(arguments.errors) +
                                                                         "_tau" + std::to_string(arguments.tau));
    if (!std::filesystem::exists(filename))
        return false;

    std::ifstream is{filename, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    iarchive(vec);
    return true;
}

} // namespace raptor

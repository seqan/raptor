// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides raptor::index_upgrader.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <raptor/argument_parsing/upgrade_arguments.hpp>
#include <raptor/index.hpp>

namespace raptor
{

class index_upgrader
{
public:
    std::string index_file{};
    std::string output_file{};
    double fpr{};
    size_t max_count{};

    index_upgrader() = default;
    index_upgrader(index_upgrader const &) = default;
    index_upgrader(index_upgrader &&) = default; // GCOVR_EXCL_LINE
    index_upgrader & operator=(index_upgrader const &) = default;
    index_upgrader & operator=(index_upgrader &&) = default;
    ~index_upgrader() = default;

    explicit index_upgrader(upgrade_arguments const & arguments, size_t const max_count) :
        index_file{arguments.index_file},
        output_file{arguments.output_file},
        fpr{arguments.fpr},
        max_count{max_count}
    {}

    void upgrade()
    {
        raptor_index<index_structure::ibf> index{};
        {
            std::ifstream is{index_file, std::ios::binary};
            cereal::BinaryInputArchive iarchive{is};
            index.load_old_index(iarchive);
        }
        if (std::isnan(fpr))
            fpr = compute_fpr(index.ibf().hash_function_count(), max_count, index.ibf().bin_size());
        index.fpr_ = fpr;
        std::cout << "FPR for " << index_file << ": " << fpr << '\n';
        index.is_hibf_ = false;
        std::ofstream os{output_file, std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(index);
    }

    static double compute_fpr(size_t const hash_fun, size_t const count, size_t const bin_size)
    {
        double const exp_arg = (hash_fun * count) / static_cast<double>(bin_size);
        double const log_arg = 1.0 - std::exp(-exp_arg);
        return std::exp(hash_fun * std::log(log_arg));
    }
};

} // namespace raptor

// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::upgrade_arguments.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <filesystem>
#include <vector>

#include <seqan3/search/kmer_index/shape.hpp>

namespace raptor
{

struct upgrade_arguments
{
    uint32_t window_size{};
    seqan3::shape shape{};
    bool compressed{};
    bool input_is_minimiser{};
    uint8_t parts{1u};
    uint8_t threads{1u};
    double fpr{std::numeric_limits<double>::quiet_NaN()};

    std::filesystem::path bin_file{};
    std::filesystem::path index_file{};
    std::filesystem::path output_file{};

    std::vector<std::vector<std::string>> bin_path{};
};

} // namespace raptor

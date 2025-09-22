// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::update_arguments.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <cstddef>    // for size_t
#include <cstdint>    // for uint8_t, uint32_t
#include <filesystem> // for path
#include <string>     // for basic_string, string
#include <vector>     // for vector

#include <seqan3/search/kmer_index/shape.hpp> // for shape, ungapped

namespace raptor
{

struct update_arguments
{
    std::filesystem::path index_file{};
    std::filesystem::path out_path{};

    std::vector<std::vector<std::string>> bin_path{};
    std::filesystem::path bin_file{};
    uint8_t threads{1u};
    bool input_is_minimiser{false};

    uint32_t window_size{20u};
    seqan3::shape shape{seqan3::ungapped{20u}};
    uint8_t shape_size{shape.size()};
    uint8_t shape_weight{shape.count()};
    uint8_t parts{1u};
    bool is_hibf{false};
    double fpr{0.05};

    std::vector<size_t> user_bins_to_delete{};
    std::vector<std::vector<std::string>> user_bins_to_insert{};
};

} // namespace raptor

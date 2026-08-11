// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::update_arguments.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <filesystem>
#include <vector>

#include <seqan3/search/kmer_index/shape.hpp>

namespace raptor
{

struct update_arguments
{
    //!\brief The index to update. Left untouched; the result is written to `out_path`.
    std::filesystem::path index_file{};
    std::filesystem::path out_path{};

    std::vector<std::vector<std::string>> bin_path{};
    //!\brief The `--insert` argument: a single sequence file, or a file containing file names.
    std::filesystem::path bin_file{};
    uint8_t threads{1u};

    // Read from the index, not from the command line.
    uint32_t window_size{20u};
    seqan3::shape shape{seqan3::ungapped{20u}};
    uint8_t shape_size{shape.size()};
    uint8_t shape_weight{shape.count()};
    uint8_t parts{1u};
    bool is_hibf{false};
    double fpr{0.05};

    /*!\brief Whether to replace every partial rebuild by a full one.
     * \details A full rebuild yields a better layout, but has to process all user bins instead of only those of one
     *          subtree.
     */
    bool disable_partial_rebuild{false};

    //!\brief The IDs of the user bins to delete.
    std::vector<size_t> user_bins_to_delete{};
    //!\brief The user bins to insert, one entry per user bin.
    std::vector<std::vector<std::string>> user_bins_to_insert{};
};

} // namespace raptor

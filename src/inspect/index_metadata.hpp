// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::detail::index_metadata and raptor::detail::read_index_metadata.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <cstddef>
#include <filesystem>
#include <optional>
#include <vector>

#include <raptor/argument_parsing/inspect_arguments.hpp>
#include <raptor/index.hpp>

namespace raptor::detail
{

//!\brief The structure of a single IBF.
struct ibf_metadata
{
    //!\brief The position of the IBF in the HIBF, or the part of an IBF index.
    size_t id{};
    //!\brief The depth of the IBF in the HIBF. `0` for the top-level IBF and for every part of an IBF index.
    size_t level{};
    //!\brief The IBF and the merged bin this IBF is the child of. Not set for a top-level IBF.
    std::optional<seqan::hibf::hierarchical_interleaved_bloom_filter::previous_ibf_id_pair> parent{};
    //!\brief The number of technical bins in use, i.e. `bin_count()`.
    size_t technical_bins{};
    //!\brief The number of technical bins the IBF can hold without being resized.
    size_t capacity{};
    //!\brief The number of merged bins, i.e. bins pointing to a lower-level IBF.
    size_t merged_bins{};
    //!\brief The number of technical bins freed by deleting a user bin.
    size_t deleted_bins{};
    //!\brief The size of a single technical bin in bits.
    size_t bin_size{};
    size_t hash_function_count{};
    //!\brief The size of the IBF in bits.
    size_t bit_size{};
};

//!\brief The structure of an index. Requires loading the whole index.
struct structure_metadata
{
    //!\brief Every IBF that is not retired, in ascending order of raptor::detail::ibf_metadata::id.
    std::vector<ibf_metadata> ibfs{};
    //!\brief The IBFs retired by an update. Their slots are kept, but they hold no data. Always empty for an IBF.
    std::vector<size_t> retired_ibfs{};
    //!\brief The user bins no technical bin refers to anymore. Always empty for an IBF.
    std::vector<size_t> deleted_user_bins{};
    //!\brief The number of levels, i.e. one more than the largest raptor::detail::ibf_metadata::level.
    size_t levels{};
};

//!\brief The metadata `raptor inspect` reports.
struct index_metadata
{
    //!\brief The path as passed to `--index`.
    std::filesystem::path index_file{};
    //!\brief The combined size of all files of the index.
    size_t size_in_bytes{};
    /*!\brief The parameters of the index, read by raptor::raptor_index::load_parameters.
     * \details
     * Its data is not loaded. raptor::raptor_index::load_parameters reads the same members for an IBF and an HIBF,
     * hence this also holds an IBF index. For an IBF index, raptor::raptor_index::config is default-constructed and
     * does not describe the index.
     */
    raptor_index<index_structure::hibf> parameters{};
    //!\brief Only set with `--detailed`.
    std::optional<structure_metadata> structure{};
};

/*!\brief Reads the metadata of the index named by `arguments`.
 * \details
 * Reads only the parameters, unless `arguments.detailed` is set. Then the whole index is loaded, one part at a time.
 * \throws sharg::validation_error if a part of a partitioned index is missing.
 * \throws sharg::parser_error if the index cannot be read or is inconsistent.
 */
index_metadata read_index_metadata(inspect_arguments const & arguments);

} // namespace raptor::detail

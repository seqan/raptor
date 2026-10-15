// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::detail::read_index_metadata.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include "index_metadata.hpp"

#include <algorithm>
#include <fstream>
#include <string>

#include <sharg/exceptions.hpp>
#include <sharg/validators.hpp>

namespace raptor::detail
{

namespace
{

//!\brief The files of the index: `<index_file>_<part>` for every part of a partitioned index, else `index_file`.
std::vector<std::filesystem::path> index_files(inspect_arguments const & arguments, size_t const parts)
{
    std::vector<std::filesystem::path> files{};

    if (arguments.is_partitioned)
    {
        files.reserve(parts);
        for (size_t part = 0; part < parts; ++part)
            files.emplace_back(arguments.index_file.string() + '_' + std::to_string(part));
    }
    else
    {
        files.push_back(arguments.index_file);
    }

    return files;
}

template <typename index_t>
index_t load_whole_index(std::filesystem::path const & path)
{
    index_t index{};
    std::ifstream stream{path, std::ios::binary};
    cereal::BinaryInputArchive archive{stream};
    archive(index);
    return index;
}

//!\brief The members of raptor::detail::ibf_metadata that do not depend on the position of `ibf` in an HIBF.
ibf_metadata describe_ibf(seqan::hibf::interleaved_bloom_filter const & ibf, size_t const id)
{
    return {.id = id,
            .technical_bins = ibf.bin_count(),
            // The IBF allocates `bin_size()` bits for every technical bin it can hold. An IBF without data, e.g. a
            // default-constructed one, has a bin size of 0.
            .capacity = ibf.bin_size() == 0u ? 0u : ibf.bit_size() / ibf.bin_size(),
            .bin_size = ibf.bin_size(),
            .hash_function_count = ibf.hash_function_count(),
            .bit_size = ibf.bit_size()};
}

//!\brief The depth of the IBF at `ibf_idx`. Follows the parents up to the top-level IBF, which is always IBF 0.
size_t level_of(seqan::hibf::hierarchical_interleaved_bloom_filter const & hibf, size_t const ibf_idx)
{
    size_t const number_of_ibfs = hibf.prev_ibf_id.size();
    size_t level{};

    for (size_t current = ibf_idx; current != 0u; ++level)
    {
        current = hibf.prev_ibf_id[current].ibf_idx;

        // A parent outside the HIBF, e.g. the marker of a retired IBF, or a cycle of parents.
        if (current >= number_of_ibfs || level >= number_of_ibfs)
            throw sharg::parser_error{"Cannot read index: IBF " + std::to_string(ibf_idx)
                                      + " is not connected to the top-level IBF."};
    }

    return level;
}

structure_metadata ibf_structure(std::vector<std::filesystem::path> const & files)
{
    structure_metadata structure{.levels = 1u};
    structure.ibfs.reserve(files.size());

    // Only one part is held in memory at a time.
    for (size_t part = 0; part < files.size(); ++part)
    {
        auto const index = load_whole_index<raptor_index<index_structure::ibf>>(files[part]);
        structure.ibfs.push_back(describe_ibf(index.ibf(), part));
    }

    return structure;
}

structure_metadata hibf_structure(std::filesystem::path const & file)
{
    auto const index = load_whole_index<raptor_index<index_structure::hibf>>(file);
    auto const & hibf = index.ibf();
    size_t const number_of_ibfs = hibf.ibf_vector.size();
    size_t const number_of_user_bins = index.bin_path().size();

    if (hibf.prev_ibf_id.size() != number_of_ibfs || hibf.ibf_bin_to_user_bin_id.size() != number_of_ibfs)
        throw sharg::parser_error{"Cannot read index: The HIBF does not store the same number of entries for every "
                                  "IBF."};

    structure_metadata structure{};
    // Deleting a user bin does not renumber the others; it only removes every reference to it.
    std::vector<bool> is_referenced(number_of_user_bins, false);

    for (size_t ibf_idx = 0; ibf_idx < number_of_ibfs; ++ibf_idx)
    {
        auto const parent = hibf.prev_ibf_id[ibf_idx];

        // See raptor::detail::retire_ibf.
        if (parent.ibf_idx == seqan::hibf::bin_kind::deleted)
        {
            structure.retired_ibfs.push_back(ibf_idx);
            continue;
        }

        ibf_metadata metadata = describe_ibf(hibf.ibf_vector[ibf_idx], ibf_idx);
        metadata.level = level_of(hibf, ibf_idx);
        if (ibf_idx != 0u)
            metadata.parent = parent;

        for (uint64_t const user_bin_id : hibf.ibf_bin_to_user_bin_id[ibf_idx])
        {
            switch (user_bin_id)
            {
            case seqan::hibf::bin_kind::merged:
                ++metadata.merged_bins;
                break;
            case seqan::hibf::bin_kind::deleted:
                ++metadata.deleted_bins;
                break;
            default:
                if (user_bin_id >= number_of_user_bins)
                    throw sharg::parser_error{"Cannot read index: IBF " + std::to_string(ibf_idx)
                                              + " refers to user bin " + std::to_string(user_bin_id)
                                              + ", but the index has only " + std::to_string(number_of_user_bins)
                                              + " user bins."};
                is_referenced[user_bin_id] = true;
            }
        }

        structure.levels = std::max(structure.levels, metadata.level + 1u);
        structure.ibfs.push_back(metadata);
    }

    for (size_t user_bin_id = 0; user_bin_id < number_of_user_bins; ++user_bin_id)
        if (!is_referenced[user_bin_id])
            structure.deleted_user_bins.push_back(user_bin_id);

    return structure;
}

} // namespace

index_metadata read_index_metadata(inspect_arguments const & arguments)
{
    std::filesystem::path const first_file =
        arguments.is_partitioned ? std::filesystem::path{arguments.index_file.string() + "_0"} : arguments.index_file;

    index_metadata metadata{.index_file = arguments.index_file};
    {
        std::ifstream stream{first_file, std::ios::binary};
        cereal::BinaryInputArchive archive{stream};
        metadata.parameters.load_parameters(archive);
    }

    auto const & parameters = metadata.parameters;

    // E.g. `--index raptor.index_0`: The parameters of a part are read fine, but the other parts would not be found.
    if (!arguments.is_partitioned && parameters.parts() > 1u)
        throw sharg::parser_error{sharg::detail::to_string("The index ",
                                                           arguments.index_file,
                                                           " is a part of a partitioned index with ",
                                                           static_cast<uint16_t>(parameters.parts()),
                                                           " parts. Pass the path without the suffix _<part>.")};

    std::vector<std::filesystem::path> const files = index_files(arguments, parameters.parts());
    sharg::input_file_validator const validator{};
    for (std::filesystem::path const & file : files)
    {
        validator(file);
        metadata.size_in_bytes += std::filesystem::file_size(file);
    }

    if (arguments.detailed)
        metadata.structure = parameters.is_hibf() ? hibf_structure(first_file) : ibf_structure(files);

    return metadata;
}

} // namespace raptor::detail

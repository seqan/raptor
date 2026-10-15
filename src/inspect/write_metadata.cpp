// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::detail::write_text and raptor::detail::write_json.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include "write_metadata.hpp"

#include <ostream>
#include <string>
#include <string_view>

#include <raptor/inspect/json_writer.hpp>

namespace raptor::detail
{

namespace
{

std::string_view index_type(raptor_index<index_structure::hibf> const & index)
{
    return index.is_hibf() ? "HIBF" : "IBF";
}

std::string_view to_string(bool const value)
{
    return value ? "true" : "false";
}

//!\brief E.g. `1,5,7`.
std::string comma_separated(std::vector<size_t> const & values)
{
    std::string result{};
    for (size_t const value : values)
    {
        if (!result.empty())
            result += ',';
        result += std::to_string(value);
    }
    return result;
}

size_t total_technical_bins(structure_metadata const & structure)
{
    size_t result{};
    for (ibf_metadata const & ibf : structure.ibfs)
        result += ibf.technical_bins;
    return result;
}

size_t total_bits(structure_metadata const & structure)
{
    size_t result{};
    for (ibf_metadata const & ibf : structure.ibfs)
        result += ibf.bit_size;
    return result;
}

} // namespace

void write_text(std::ostream & stream, index_metadata const & metadata, bool const user_bins)
{
    auto const & index = metadata.parameters;
    auto const & config = index.config();
    seqan3::shape const shape = index.shape();

    stream << "### Index\n";
    // Quoted, as in the header of a `raptor search` result.
    stream << "## Index = " << metadata.index_file << '\n';
    stream << "## Index version = " << raptor_index<>::version << '\n';
    stream << "## Index type = " << index_type(index) << '\n';
    stream << "## Index parts = " << static_cast<uint16_t>(index.parts()) << '\n';
    stream << "## Index size [bytes] = " << metadata.size_in_bytes << '\n';

    stream << "### Minimiser parameters\n";
    stream << "## Window size = " << index.window_size() << '\n';
    stream << "## Shape = " << shape.to_string() << '\n';
    // The size type of seqan3::shape is a `uint8_t`, which would be printed as a character.
    stream << "## Shape size (length) = " << static_cast<size_t>(shape.size()) << '\n';
    stream << "## Shape count (number of 1s) = " << static_cast<size_t>(shape.count()) << '\n';

    stream << "### Bloom filter parameters\n";
    stream << "## False positive rate = " << shortest_double(index.fpr()) << '\n';

    if (index.is_hibf())
    {
        stream << "## Relaxed false positive rate = " << shortest_double(config.relaxed_fpr) << '\n';
        stream << "## Hash functions = " << config.number_of_hash_functions << '\n';

        stream << "### Layout parameters\n";
        stream << "## Number of user bins in layout = " << config.number_of_user_bins << '\n';
        stream << "## tmax = " << config.tmax << '\n';
        stream << "## Empty bin fraction = " << shortest_double(config.empty_bin_fraction) << '\n';
        stream << "## Sketch bits = " << static_cast<uint16_t>(config.sketch_bits) << '\n';
        stream << "## alpha = " << shortest_double(config.alpha) << '\n';
        stream << "## Max rearrangement ratio = " << shortest_double(config.max_rearrangement_ratio) << '\n';
        stream << "## Disable estimate union = " << to_string(config.disable_estimate_union) << '\n';
        stream << "## Disable rearrangement = " << to_string(config.disable_rearrangement) << '\n';

        stream << "### Update state\n";
        stream << "## Track occupancy = " << to_string(config.track_occupancy) << '\n';
        stream << "## Number of IBFs = " << index.number_of_ibfs() << '\n';
        stream << "## Original number of IBFs = " << index.original_number_of_ibfs() << '\n';
        stream << "## Resized IBFs = " << index.number_of_resized_ibfs() << '\n';
    }

    stream << "### User bins\n";
    stream << "## Number of user bins = " << index.bin_path().size() << '\n';

    // The same format as the user bins in the header of a `raptor search` result.
    if (user_bins)
    {
        size_t user_bin_id{};
        for (std::vector<std::string> const & files : index.bin_path())
        {
            stream << '#' << user_bin_id << '\t';
            char separator{};
            for (std::string const & file : files)
            {
                if (separator)
                    stream << separator;
                stream << file;
                separator = ',';
            }
            stream << '\n';
            ++user_bin_id;
        }
    }

    if (metadata.structure)
    {
        structure_metadata const & structure = *metadata.structure;
        stream << "### Structure\n";
        stream << "## Levels = " << structure.levels << '\n';
        stream << "## Technical bins = " << total_technical_bins(structure) << '\n';
        stream << "## Size [bits] = " << total_bits(structure) << '\n';

        if (index.is_hibf())
        {
            stream << "## Number of retired IBFs = " << structure.retired_ibfs.size() << '\n';
            if (!structure.retired_ibfs.empty())
                stream << "## Retired IBFs = " << comma_separated(structure.retired_ibfs) << '\n';
            stream << "## Number of deleted user bins = " << structure.deleted_user_bins.size() << '\n';
            if (!structure.deleted_user_bins.empty())
                stream << "## Deleted user bins = " << comma_separated(structure.deleted_user_bins) << '\n';
        }

        stream << "#IBF\tLEVEL\tPARENT_IBF\tPARENT_BIN\tTECHNICAL_BINS\tCAPACITY\tMERGED_BINS\tDELETED_BINS\tBIN_SIZE"
                  "\tHASH_FUNCTIONS\tBITS\n";

        for (ibf_metadata const & ibf : structure.ibfs)
        {
            stream << ibf.id << '\t' << ibf.level << '\t';
            if (ibf.parent)
                stream << ibf.parent->ibf_idx << '\t' << ibf.parent->bin_idx << '\t';
            else
                stream << "-\t-\t";
            stream << ibf.technical_bins << '\t' << ibf.capacity << '\t' << ibf.merged_bins << '\t' << ibf.deleted_bins
                   << '\t' << ibf.bin_size << '\t' << ibf.hash_function_count << '\t' << ibf.bit_size << '\n';
        }
    }
}

void write_json(std::ostream & stream, index_metadata const & metadata, bool const user_bins)
{
    auto const & index = metadata.parameters;
    auto const & config = index.config();
    seqan3::shape const shape = index.shape();
    json_writer json{stream};

    json.begin_object();

    json.key("index").begin_object();
    json.key("path").string(metadata.index_file.string());
    json.key("version").number(raptor_index<>::version);
    json.key("type").string(index_type(index));
    json.key("parts").number(index.parts());
    json.key("size_in_bytes").number(metadata.size_in_bytes);
    json.end_object();

    json.key("minimiser").begin_object();
    json.key("window_size").number(index.window_size());
    json.key("shape").string(shape.to_string());
    json.key("shape_size").number(shape.size());
    json.key("shape_count").number(shape.count());
    json.end_object();

    json.key("bloom_filter").begin_object();
    json.key("false_positive_rate").number(index.fpr());
    if (index.is_hibf())
    {
        json.key("relaxed_false_positive_rate").number(config.relaxed_fpr);
        json.key("hash_functions").number(config.number_of_hash_functions);
    }
    json.end_object();

    if (index.is_hibf())
    {
        json.key("layout").begin_object();
        json.key("number_of_user_bins").number(config.number_of_user_bins);
        json.key("tmax").number(config.tmax);
        json.key("empty_bin_fraction").number(config.empty_bin_fraction);
        json.key("sketch_bits").number(config.sketch_bits);
        json.key("alpha").number(config.alpha);
        json.key("max_rearrangement_ratio").number(config.max_rearrangement_ratio);
        json.key("disable_estimate_union").boolean(config.disable_estimate_union);
        json.key("disable_rearrangement").boolean(config.disable_rearrangement);
        json.end_object();

        json.key("update").begin_object();
        json.key("track_occupancy").boolean(config.track_occupancy);
        json.key("number_of_ibfs").number(index.number_of_ibfs());
        json.key("original_number_of_ibfs").number(index.original_number_of_ibfs());
        json.key("resized_ibfs").number(index.number_of_resized_ibfs());
        json.end_object();
    }

    json.key("user_bins").begin_object();
    json.key("count").number(index.bin_path().size());
    if (user_bins)
    {
        // The position in the array is the user bin ID.
        json.key("paths").begin_array();
        for (std::vector<std::string> const & files : index.bin_path())
        {
            json.begin_array(true);
            for (std::string const & file : files)
                json.string(file);
            json.end_array();
        }
        json.end_array();
    }
    json.end_object();

    if (metadata.structure)
    {
        structure_metadata const & structure = *metadata.structure;

        json.key("structure").begin_object();
        json.key("levels").number(structure.levels);
        json.key("technical_bins").number(total_technical_bins(structure));
        json.key("bits").number(total_bits(structure));

        if (index.is_hibf())
        {
            json.key("retired_ibfs").begin_array(true);
            for (size_t const ibf_idx : structure.retired_ibfs)
                json.number(ibf_idx);
            json.end_array();

            json.key("deleted_user_bins").begin_array(true);
            for (size_t const user_bin_id : structure.deleted_user_bins)
                json.number(user_bin_id);
            json.end_array();
        }

        json.key("ibfs").begin_array();
        for (ibf_metadata const & ibf : structure.ibfs)
        {
            json.begin_object(true);
            json.key("id").number(ibf.id);
            json.key("level").number(ibf.level);
            if (ibf.parent)
            {
                json.key("parent_ibf").number(ibf.parent->ibf_idx);
                json.key("parent_bin").number(ibf.parent->bin_idx);
            }
            else
            {
                json.key("parent_ibf").null();
                json.key("parent_bin").null();
            }
            json.key("technical_bins").number(ibf.technical_bins);
            json.key("capacity").number(ibf.capacity);
            json.key("merged_bins").number(ibf.merged_bins);
            json.key("deleted_bins").number(ibf.deleted_bins);
            json.key("bin_size").number(ibf.bin_size);
            json.key("hash_functions").number(ibf.hash_function_count);
            json.key("bits").number(ibf.bit_size);
            json.end_object();
        }
        json.end_array();

        json.end_object();
    }

    json.end_object();
}

} // namespace raptor::detail

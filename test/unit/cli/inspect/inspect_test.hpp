// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides a fixture for `raptor inspect` CLI tests.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <cereal/archives/json.hpp> // Configures the rapidjson bundled with cereal.

#include "../update/update_test.hpp"

/*!\brief Fixture for `raptor inspect` tests.
 * \details
 * The output of `raptor inspect` is checked against the index loaded via the API, which does not share code with
 * the implementation of `raptor inspect`. Derives from update_test to build indices that can be updated.
 */
struct inspect_test : public update_test
{
    using ibf_index_t = raptor::raptor_index<raptor::index_structure::ibf>;
    using hibf_index_t = raptor::raptor_index<raptor::index_structure::hibf>;

    //!\brief The text output of `raptor inspect`, split into its parts.
    struct text_output
    {
        //!\brief The `### Section` lines, without `### `.
        std::vector<std::string> sections{};
        //!\brief The `## Key = value` lines.
        std::map<std::string, std::string> values{};
        //!\brief The `#<id>\t<files>` lines, in order.
        std::vector<std::string> user_bins{};
        //!\brief The `#IBF\t...` line.
        std::string ibf_columns{};
        //!\brief The tab-separated lines after `#IBF\t...`.
        std::vector<std::vector<std::string>> ibfs{};
    };

    //!\brief An IBF as `raptor inspect --detailed` should report it, computed from the index loaded via the API.
    struct expected_ibf
    {
        size_t id{};
        size_t level{};
        std::string parent_ibf{"-"};
        std::string parent_bin{"-"};
        size_t technical_bins{};
        size_t capacity{};
        size_t merged_bins{};
        size_t deleted_bins{};
        size_t bin_size{};
        size_t hash_function_count{};
        size_t bit_size{};

        std::vector<std::string> as_row() const
        {
            return {std::to_string(id),
                    std::to_string(level),
                    parent_ibf,
                    parent_bin,
                    std::to_string(technical_bins),
                    std::to_string(capacity),
                    std::to_string(merged_bins),
                    std::to_string(deleted_bins),
                    std::to_string(bin_size),
                    std::to_string(hash_function_count),
                    std::to_string(bit_size)};
        }
    };

    static inline std::string const ibf_columns{"#IBF\tLEVEL\tPARENT_IBF\tPARENT_BIN\tTECHNICAL_BINS\tCAPACITY\t"
                                                "MERGED_BINS\tDELETED_BINS\tBIN_SIZE\tHASH_FUNCTIONS\tBITS"};

    static std::vector<std::string> split(std::string const & line, char const delimiter)
    {
        std::vector<std::string> result{};
        size_t start{};
        for (size_t end = line.find(delimiter); end != std::string::npos; end = line.find(delimiter, start))
        {
            result.push_back(line.substr(start, end - start));
            start = end + 1u;
        }
        result.push_back(line.substr(start));
        return result;
    }

    static text_output parse_text(std::string const & text)
    {
        text_output result{};
        std::istringstream stream{text};
        bool in_table{false};

        for (std::string line{}; std::getline(stream, line);)
        {
            if (in_table)
            {
                result.ibfs.push_back(split(line, '\t'));
            }
            else if (line.starts_with("### "))
            {
                result.sections.push_back(line.substr(4u));
            }
            else if (line.starts_with("## "))
            {
                size_t const separator = line.find(" = ");
                EXPECT_NE(separator, std::string::npos) << line;
                std::string const key = line.substr(3u, separator - 3u);
                EXPECT_FALSE(result.values.contains(key)) << "Duplicate key: " << key;
                result.values[key] = line.substr(separator + 3u);
            }
            else if (line.starts_with("#IBF\t"))
            {
                result.ibf_columns = line;
                in_table = true;
            }
            else
            {
                EXPECT_TRUE(line.starts_with('#')) << line;
                result.user_bins.push_back(line);
            }
        }

        return result;
    }

    //!\brief The user bin lines as in the header of a `raptor search` result.
    static std::vector<std::string> user_bin_lines(std::vector<std::vector<std::string>> const & bin_path)
    {
        std::vector<std::string> result{};
        for (size_t user_bin_id = 0; user_bin_id < bin_path.size(); ++user_bin_id)
        {
            std::string line{'#' + std::to_string(user_bin_id) + '\t'};
            for (size_t i = 0; i < bin_path[user_bin_id].size(); ++i)
                line += (i ? "," : "") + bin_path[user_bin_id][i];
            result.push_back(line);
        }
        return result;
    }

    // Strict RFC 8259: Unlike cereal's default flags, this rejects NaN and infinity.
    static rapidjson::Document parse_json(std::string const & text)
    {
        rapidjson::Document document{};
        document.Parse<rapidjson::kParseFullPrecisionFlag | rapidjson::kParseValidateEncodingFlag>(text.data(),
                                                                                                   text.size());
        EXPECT_FALSE(document.HasParseError()) << "Invalid JSON at offset " << document.GetErrorOffset() << ":\n"
                                               << text;
        return document;
    }

    template <typename index_t>
    static index_t load_index(std::filesystem::path const & path)
    {
        index_t index{};
        std::ifstream stream{path, std::ios::binary};
        cereal::BinaryInputArchive archive{stream};
        archive(index);
        return index;
    }

    template <typename index_t>
    static void store_index(std::filesystem::path const & path, index_t const & index)
    {
        std::ofstream stream{path, std::ios::binary};
        cereal::BinaryOutputArchive archive{stream};
        archive(index);
    }

    static expected_ibf expected_from(seqan::hibf::interleaved_bloom_filter const & ibf, size_t const id)
    {
        return {.id = id,
                .technical_bins = ibf.bin_count(),
                .capacity = ibf.occupancy.size(),
                .bin_size = ibf.bin_size(),
                .hash_function_count = ibf.hash_function_count(),
                .bit_size = ibf.bit_size()};
    }

    //!\brief Whether `prev_ibf_id` marks the IBF as retired, see raptor::detail::retire_ibf.
    static bool is_retired(hibf_index_t const & index, size_t const ibf_idx)
    {
        return index.ibf().prev_ibf_id[ibf_idx].ibf_idx == seqan::hibf::bin_kind::deleted;
    }

    //!\brief Every IBF that is not retired.
    static std::vector<expected_ibf> expected_ibfs(hibf_index_t const & index)
    {
        auto const & hibf = index.ibf();
        std::vector<expected_ibf> result{};

        for (size_t ibf_idx = 0; ibf_idx < hibf.ibf_vector.size(); ++ibf_idx)
        {
            if (is_retired(index, ibf_idx))
                continue;

            // `occupancy` always has one entry per technical bin the IBF can hold.
            expected_ibf ibf = expected_from(hibf.ibf_vector[ibf_idx], ibf_idx);

            for (size_t current = ibf_idx; current != 0u; current = hibf.prev_ibf_id[current].ibf_idx)
                ++ibf.level;

            if (ibf_idx != 0u)
            {
                ibf.parent_ibf = std::to_string(hibf.prev_ibf_id[ibf_idx].ibf_idx);
                ibf.parent_bin = std::to_string(hibf.prev_ibf_id[ibf_idx].bin_idx);
            }

            ibf.merged_bins = std::ranges::count(hibf.ibf_bin_to_user_bin_id[ibf_idx], seqan::hibf::bin_kind::merged);
            ibf.deleted_bins = std::ranges::count(hibf.ibf_bin_to_user_bin_id[ibf_idx], seqan::hibf::bin_kind::deleted);
            result.push_back(ibf);
        }

        return result;
    }

    static std::vector<size_t> expected_retired_ibfs(hibf_index_t const & index)
    {
        std::vector<size_t> result{};
        for (size_t ibf_idx = 0; ibf_idx < index.ibf().ibf_vector.size(); ++ibf_idx)
            if (is_retired(index, ibf_idx))
                result.push_back(ibf_idx);
        return result;
    }

    //!\brief The user bins that no technical bin of an IBF that is not retired refers to.
    static std::vector<size_t> expected_deleted_user_bins(hibf_index_t const & index)
    {
        std::set<size_t> referenced{};
        for (size_t ibf_idx = 0; ibf_idx < index.ibf().ibf_vector.size(); ++ibf_idx)
            if (!is_retired(index, ibf_idx))
                for (uint64_t const user_bin_id : index.ibf().ibf_bin_to_user_bin_id[ibf_idx])
                    if (user_bin_id < index.bin_path().size())
                        referenced.insert(user_bin_id);

        std::vector<size_t> result{};
        for (size_t user_bin_id = 0; user_bin_id < index.bin_path().size(); ++user_bin_id)
            if (!referenced.contains(user_bin_id))
                result.push_back(user_bin_id);
        return result;
    }

    static std::string comma_separated(std::vector<size_t> const & values)
    {
        std::string result{};
        for (size_t const value : values)
            result += (result.empty() ? "" : ",") + std::to_string(value);
        return result;
    }

    static std::vector<size_t> json_ids(rapidjson::Value const & array)
    {
        std::vector<size_t> result{};
        for (auto const & value : array.GetArray())
            result.push_back(value.GetUint64());
        return result;
    }

    //!\brief Checks the `ibfs` array of the JSON output against `expected`.
    static void expect_json_ibfs(rapidjson::Value const & ibfs, std::vector<expected_ibf> const & expected)
    {
        ASSERT_EQ(ibfs.Size(), expected.size());
        for (size_t i = 0; i < expected.size(); ++i)
        {
            rapidjson::Value const & ibf = ibfs[static_cast<rapidjson::SizeType>(i)];
            auto const parent = [&](char const * name) -> std::string
            {
                return ibf[name].IsNull() ? "-" : std::to_string(ibf[name].GetUint64());
            };

            EXPECT_EQ(ibf["id"].GetUint64(), expected[i].id);
            EXPECT_EQ(ibf["level"].GetUint64(), expected[i].level);
            EXPECT_EQ(parent("parent_ibf"), expected[i].parent_ibf);
            EXPECT_EQ(parent("parent_bin"), expected[i].parent_bin);
            EXPECT_EQ(ibf["technical_bins"].GetUint64(), expected[i].technical_bins);
            EXPECT_EQ(ibf["capacity"].GetUint64(), expected[i].capacity);
            EXPECT_EQ(ibf["merged_bins"].GetUint64(), expected[i].merged_bins);
            EXPECT_EQ(ibf["deleted_bins"].GetUint64(), expected[i].deleted_bins);
            EXPECT_EQ(ibf["bin_size"].GetUint64(), expected[i].bin_size);
            EXPECT_EQ(ibf["hash_functions"].GetUint64(), expected[i].hash_function_count);
            EXPECT_EQ(ibf["bits"].GetUint64(), expected[i].bit_size);
        }
    }

    //!\brief Checks the rows of the text output against `expected`.
    static void expect_text_ibfs(text_output const & output, std::vector<expected_ibf> const & expected)
    {
        EXPECT_EQ(output.ibf_columns, ibf_columns);
        ASSERT_EQ(output.ibfs.size(), expected.size());
        for (size_t i = 0; i < expected.size(); ++i)
            EXPECT_EQ(output.ibfs[i], expected[i].as_row()) << "Row " << i;
    }
};

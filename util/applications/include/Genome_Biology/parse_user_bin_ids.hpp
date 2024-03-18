// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <charconv>
#include <filesystem>
#include <fstream>
#include <string>
#include <string_view>

#include <hibf/contrib/robin_hood.hpp>

inline robin_hood::unordered_map<std::string, uint64_t>
parse_user_bin_ids(std::filesystem::path const & user_bin_ids_file)
{
    std::string line_buffer{};
    uint64_t id_buffer{};
    robin_hood::unordered_map<std::string, uint64_t> ub_name_to_id;
    std::ifstream user_bin_ids_in{user_bin_ids_file};

    // Contains lines: "some_number <tab> reference_name"
    while (std::getline(user_bin_ids_in, line_buffer))
    {
        auto tab_it{line_buffer.begin() + line_buffer.find('\t')};
        std::string_view const id_value{line_buffer.begin(), tab_it};
        std::string_view const name_key{++tab_it, line_buffer.end()};
        std::from_chars(id_value.data(), id_value.data() + id_value.size(), id_buffer);
        ub_name_to_id.emplace(name_key, id_buffer);
    }
    return ub_name_to_id;
}

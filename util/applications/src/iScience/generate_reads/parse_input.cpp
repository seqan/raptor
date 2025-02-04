// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <charconv>
#include <filesystem>
#include <fstream>
#include <string>
#include <tuple>

#include <sharg/exceptions.hpp>
#include <sharg/validators.hpp>

#include "cmd_arguments.hpp"

namespace raptor::util::generate_reads
{

void parse_input(cmd_arguments & arguments)
{
    std::ifstream istrm{arguments.bin_file_path};
    std::string line;
    sharg::input_file_validator validator{};
    std::filesystem::path bin_path{};

    while (std::getline(istrm, line))
    {
        if (line.empty())
            continue;

        auto const & [path_span, weight_span] = [&line]()
        {
            auto const tab = std::ranges::find(line, '\t');
            auto const end = line.end();

            std::span path_span{line.begin(), tab};
            std::span weight_span{std::ranges::next(tab, 1, end), end};

            return std::make_tuple(path_span, weight_span);
        }();

        bin_path.assign(path_span.begin(), path_span.end());
        validator(bin_path);
        arguments.bin_path.push_back(bin_path);
        ++arguments.number_of_bins;

        if (arguments.weight_mode != weight::from_weight_column)
            continue;

        if (!weight_span.empty())
        {
            arguments.number_of_reads_per_bin.emplace_back();
            std::from_chars(weight_span.data(),
                            weight_span.data() + weight_span.size(),
                            arguments.number_of_reads_per_bin.back());
        }
        else
        {
            throw sharg::parser_error{"No weight given for bin: " + bin_path.string()};
        }
    }

    if (arguments.weight_mode != weight::from_weight_column)
        arguments.number_of_reads_per_bin.resize(arguments.number_of_bins);
}

} // namespace raptor::util::generate_reads

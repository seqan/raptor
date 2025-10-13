// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::parse_bin_path.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <filesystem> // for path, operator==
#include <fstream>    // for basic_istream, basic_ifstream, basic_ios, getline
#include <sstream>    // for basic_stringstream
#include <string>     // for basic_string, char_traits, string
#include <utility>    // for move
#include <vector>     // for vector

#include <chopper/layout/input.hpp> // for read_filenames_from

#include <raptor/argument_parsing/build_arguments.hpp>   // for build_arguments
#include <raptor/argument_parsing/parse_bin_path.hpp>    // for parse_bin_path
#include <raptor/argument_parsing/prepare_arguments.hpp> // for prepare_arguments
#include <raptor/argument_parsing/validators.hpp>        // for bin_validator

namespace raptor
{

namespace detail
{

void parse_bin_path(std::filesystem::path const & bin_file,
                    std::vector<std::vector<std::string>> & bin_path,
                    bool const is_hibf)
{
    std::ifstream istrm{bin_file};

    if (is_hibf)
    {
        bin_path = chopper::layout::read_filenames_from(istrm);
    }
    else
    {
        std::string line{};
        std::string file_name{};
        std::vector<std::string> tmp{};

        while (std::getline(istrm, line))
        {
            if (!line.empty())
            {
                tmp.clear();
                std::stringstream sstream{line};

                while (std::getline(sstream, file_name, ' '))
                    if (!file_name.empty())
                        tmp.emplace_back(std::move(file_name));

                bin_path.emplace_back(std::move(tmp));
            }
        }
    }

    bin_validator{}(bin_path);
}

} // namespace detail

void parse_bin_path(build_arguments & arguments)
{
    raptor::detail::parse_bin_path(arguments.bin_file, arguments.bin_path, arguments.is_hibf);
    arguments.bins = arguments.bin_path.size();
    std::filesystem::path first_bin_path = arguments.bin_path[0][0];
    arguments.input_is_minimiser = first_bin_path.extension() == ".minimiser";
}

void parse_bin_path(prepare_arguments & arguments)
{
    raptor::detail::parse_bin_path(arguments.bin_file, arguments.bin_path, false);
}

} // namespace raptor

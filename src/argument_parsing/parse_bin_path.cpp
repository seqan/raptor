// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <raptor/argument_parsing/parse_bin_path.hpp>
#include <raptor/argument_parsing/validators.hpp>
#include <raptor/build/hibf/parse_chopper_pack_line.hpp>

namespace raptor
{

namespace detail
{

void parse_bin_path(std::filesystem::path const & bin_file,
                    std::vector<std::vector<std::string>> & bin_path,
                    bool const is_hibf)
{
    std::ifstream istrm{bin_file};
    std::string line{};
    std::string file_name{};
    std::vector<std::string> tmp{};

    if (is_hibf)
    {
        while (std::getline(istrm, line) && line.substr(0, 6) != "#FILES")
        {}
        while (std::getline(istrm, line))
            bin_path.push_back(std::move(hibf::parse_chopper_pack_line(line).filenames));
    }
    else
    {
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
    raptor::detail::parse_bin_path(arguments.bin_file, arguments.original_bin_path, arguments.is_hibf);
    arguments.bins = arguments.original_bin_path.size();
    std::filesystem::path first_bin_path = arguments.original_bin_path[0][0];
    arguments.input_is_minimiser = first_bin_path.extension() == ".minimiser";
}

void parse_bin_path(prepare_arguments & arguments)
{
    raptor::detail::parse_bin_path(arguments.bin_file, arguments.bin_path, false);
}

void parse_bin_path(upgrade_arguments & arguments)
{
    raptor::detail::parse_bin_path(arguments.bin_file, arguments.bin_path, false);
}

} // namespace raptor

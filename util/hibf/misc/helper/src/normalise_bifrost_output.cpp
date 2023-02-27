// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <cassert>
#include <charconv>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <robin_hood.h>

#include <raptor/argument_parsing/validators.hpp>

inline robin_hood::unordered_map<std::string, std::string>
parse_user_bin_ids(std::filesystem::path const & user_bin_ids_file)
{
    std::string line_buffer{};
    robin_hood::unordered_map<std::string, std::string> ub_name_to_id;
    std::ifstream user_bin_ids_in{user_bin_ids_file};

    // Contains lines: "some_number <tab> reference_name"
    while (std::getline(user_bin_ids_in, line_buffer))
    {
        auto tab_it{line_buffer.begin() + line_buffer.find('\t')};
        std::string_view const id_value{line_buffer.begin(), tab_it};
        std::string_view const name_key{++tab_it, line_buffer.end()};
        ub_name_to_id.emplace(name_key, id_value);
    }
    return ub_name_to_id;
}

inline void check_output_file(std::filesystem::path const & output_file)
{
    std::filesystem::path const output_directory = output_file.parent_path();
    std::error_code ec{};
    std::filesystem::create_directories(output_directory, ec);

    if (!output_directory.empty() && ec)
        throw seqan3::argument_parser_error{
            seqan3::detail::to_string("Failed to create directory\"", output_directory.c_str(), "\": ", ec.message())};
}

struct config
{
    std::filesystem::path query_names_file{};
    std::filesystem::path user_bin_ids_file{};
    std::filesystem::path bifrost_result_file{};
    std::filesystem::path output_file{};
};

std::vector<std::string> parse_query_names(std::filesystem::path const & query_names_file)
{
    std::string line_buffer{};
    std::vector<std::string> query_names;
    std::ifstream query_names_in{query_names_file};

    std::cerr << "Reading " << query_names_file << " ... " << std::flush;
    // Contains lines: "query_name"
    while (std::getline(query_names_in, line_buffer))
        query_names.push_back(line_buffer);
    std::cerr << "Done" << std::endl;
    return query_names;
}

void normalise_output(config const & cfg)
{
    // All query names
    std::vector<std::string> const query_names{parse_query_names(cfg.query_names_file)};
    // map[reference_name] = number
    std::cerr << "Reading " << cfg.user_bin_ids_file << " ... " << std::flush;
    robin_hood::unordered_map<std::string, std::string> const ub_name_to_id{parse_user_bin_ids(cfg.user_bin_ids_file)};
    std::cerr << "Done" << std::endl;

    // Process bifrost results
    std::ifstream bifrost_result_in{cfg.bifrost_result_file};
    std::ofstream bifrost_result_out{cfg.output_file};
    size_t current_query_number{};
    std::vector<uint64_t> results;

    // Buffers for file I/O
    std::string ub_name_buffer{};
    std::string result_buffer{};
    std::string line_buffer{};

    std::vector<std::string> bifrost_user_bins{};

    std::string normalised_bifrost_line{};

    // ## Bifrost results ##
    // Bifrost outputs a matrix. Column names = user bin id. Row names = query names
    auto split_line_by_tab_and = [](std::string_view bifrost_line, auto do_me)
    {
        std::string_view::size_type current_pos = 0;
        std::string_view::size_type tab_pos{bifrost_line.find('\t')};
        size_t column_idx{};

        while (tab_pos != std::string_view::npos)
        {
            auto current = std::string(&bifrost_line[current_pos], tab_pos - current_pos);
            do_me(current, column_idx);
            current_pos = tab_pos + 1;
            tab_pos = bifrost_line.find('\t', current_pos);
            ++column_idx;
        }
        // process last ub
        auto last = std::string(&bifrost_line[current_pos], bifrost_line.size() - current_pos);
        do_me(last, column_idx);
    };

    auto parse_header_user_bin_id = [&bifrost_user_bins, &ub_name_to_id](std::string const & sv, size_t idx)
    {
        if (idx != 0)
        {
            auto filename_start = sv.find_last_of('/') + 1;
            auto user_bin_id = sv.substr(filename_start, sv.size() - filename_start - 7 /* |".fna.gz"| */);

            try
            {
                bifrost_user_bins.push_back(ub_name_to_id.at(user_bin_id));
            }
            catch (std::exception const & e)
            {
                std::cerr << "Could not find id: " << user_bin_id << std::endl;
                throw e;
            }
        }
        else
        {
            bifrost_user_bins.push_back("0"); // don't mess up the indices
        }
    };

    auto insert_if_one = [&normalised_bifrost_line, &bifrost_user_bins](std::string_view sv, size_t idx)
    {
        if (sv == std::string_view{"1"}) // excludes 0 and the first column which is alywas the query name
        {
            normalised_bifrost_line.insert(normalised_bifrost_line.end(),
                                           bifrost_user_bins[idx].begin(),
                                           bifrost_user_bins[idx].end());
            normalised_bifrost_line.push_back(',');
        }
    };

    std::cerr << "Processing " << cfg.bifrost_result_file << " ... " << std::endl;

    // Parse header line
    if (std::getline(bifrost_result_in, line_buffer))
    {
        assert(line_buffer.starts_with("query_name"));
        split_line_by_tab_and(line_buffer, parse_header_user_bin_id);
    }

    std::cerr << "Successfully parsed Header line ... " << std::endl;
    assert(ub_name_to_id.size() == bifrost_user_bins.size());

    while (std::getline(bifrost_result_in, line_buffer))
    {
        assert(query_names[current_query_number] == std::string(&line_buffer[0], line_buffer.find('\t')));
        normalised_bifrost_line = query_names[current_query_number];
        normalised_bifrost_line.push_back('\t');
        split_line_by_tab_and(line_buffer, insert_if_one);
        if (normalised_bifrost_line.back() == ',')
            normalised_bifrost_line.pop_back();
        bifrost_result_out << normalised_bifrost_line << '\n';
        ++current_query_number;
        normalised_bifrost_line.clear();
    }

    std::cerr << "Done" << std::endl;
}

void init_parser(seqan3::argument_parser & parser, config & cfg)
{
    parser.add_option(cfg.query_names_file,
                      '\0',
                      "query_names",
                      "The file containing query names, e.g., \"query.names\".",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{});
    parser.add_option(cfg.user_bin_ids_file,
                      '\0',
                      "user_bin_ids",
                      "The file containing user bin ids, e.g., \"user_bin.ids\".",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{});
    parser.add_option(cfg.bifrost_result_file,
                      '\0',
                      "bifrost_results",
                      "The bifrost result file, e.g., \"bifrost.results\".",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{});
    parser.add_option(cfg.output_file,
                      '\0',
                      "output_file",
                      "Provide a path to the output.",
                      seqan3::option_spec::required);
}

int main(int argc, char ** argv)
{
    seqan3::argument_parser parser{"normalise_bifrost_output", argc, argv, seqan3::update_notifications::off};
    parser.info.author = "Svenja Mehringer, Enrico Seiler";
    parser.info.email = "enrico.seiler@fu-berlin.de";
    parser.info.short_description = "Converts bifrost results into raptor-like results.";
    parser.info.version = "0.0.1";

    config cfg{};
    init_parser(parser, cfg);

    try
    {
        parser.parse();
        check_output_file(cfg.output_file);
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        std::cerr << "[Error] " << ext.what() << '\n';
        std::exit(-1);
    }

    normalise_output(cfg);
}

// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <algorithm>
#include <cassert>
#include <charconv>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <robin_hood.h>

#include <raptor/argument_parsing/validators.hpp>

inline void check_output_file(std::filesystem::path const & output_file)
{
    std::filesystem::path const output_directory = output_file.parent_path();
    std::error_code ec{};
    std::filesystem::create_directories(output_directory, ec);

    if (!output_directory.empty() && ec)
        throw seqan3::argument_parser_error{
            seqan3::detail::to_string("Failed to create directory\"", output_directory.c_str(), "\": ", ec.message())};
}

enum class validation
{
    fps,
    fns
};

inline robin_hood::unordered_map<std::string, std::vector<uint64_t>>
parse_validation_for(std::filesystem::path const & validation_file, validation which)
{
    uint64_t user_bin_id_buffer{};
    uint64_t count_buffer{};
    robin_hood::unordered_map<std::string, std::vector<uint64_t>> truths_per_query;

    std::ifstream validation_in{validation_file};
    std::string line_buffer{};

    // Contains lines: "squery_name:user_in_id:count"
    while (std::getline(validation_in, line_buffer))
    {
        auto start{line_buffer.begin()};
        auto colon{start + line_buffer.find(':')};
        std::string const query_name{line_buffer.begin(), colon};
        start = ++colon;
        colon = line_buffer.begin() + line_buffer.find(':', colon - line_buffer.begin());
        std::string_view const user_in_id{start, colon};
        std::string_view const count_str{++colon, line_buffer.end()};

        std::from_chars(user_in_id.data(), user_in_id.data() + user_in_id.size(), user_bin_id_buffer);
        std::from_chars(count_str.data(), count_str.data() + count_str.size(), count_buffer);

        if ((which == validation::fps && count_buffer < 155) || (which == validation::fns && count_buffer >= 155))
            truths_per_query[query_name].push_back(user_bin_id_buffer);
    }

    return truths_per_query;
}

struct config
{
    std::filesystem::path raptor_result_file{};
    std::filesystem::path FN_file{};
    std::filesystem::path FP_file{};
    std::filesystem::path output_file{};
};

void correct_truth_file(config const & cfg)
{
    // Process raptor results
    std::ifstream raptor_result_in{cfg.raptor_result_file};
    std::ofstream raptor_result_out{cfg.output_file};

    // Buffers for file I/O
    std::vector<uint64_t> result_user_bins{};
    std::string line_buffer{};
    uint64_t user_bin_id_buffer{};

    auto fns_truths_per_query = parse_validation_for(cfg.FN_file, validation::fns);
    auto fps_truths_per_query = parse_validation_for(cfg.FP_file, validation::fps);

    std::cerr << "Processsssssing " << cfg.raptor_result_file << " ... " << std::flush;

    // rewrite header
    while (std::getline(raptor_result_in, line_buffer) && line_buffer[0] == '#')
        raptor_result_out << line_buffer << '\n';

    do
    {
        result_user_bins.clear();

        auto tab_it{line_buffer.begin() + line_buffer.find('\t')};
        std::string const id{line_buffer.begin(), tab_it};
        std::string_view const bins{++tab_it, line_buffer.end()};

        if (bins.empty())
        {
            raptor_result_out << id << '\t' << '\n';
            continue;
        }

        std::string_view::size_type current_pos = 0;
        std::string_view::size_type comma_pos{bins.find(',')};

        // get all user bins
        while (comma_pos != std::string_view::npos)
        {
            auto user_bin_id = std::string(&bins[current_pos], comma_pos - current_pos);
            std::from_chars(user_bin_id.data(), user_bin_id.data() + user_bin_id.size(), user_bin_id_buffer);
            result_user_bins.push_back(user_bin_id_buffer);
            current_pos = comma_pos + 1;
            comma_pos = bins.find(',', current_pos);
        }
        // process last ub
        auto user_bin_id = std::string(&bins[current_pos], bins.size() - current_pos);
        std::from_chars(user_bin_id.data(), user_bin_id.data() + user_bin_id.size(), user_bin_id_buffer);
        result_user_bins.push_back(user_bin_id_buffer);

        // remove true false posiives
        auto find_fps = [&](auto s)
        {
            auto & l = fps_truths_per_query[id];
            return std::find(l.begin(), l.end(), s) != l.end();
        };
        result_user_bins.erase(std::remove_if(result_user_bins.begin(), result_user_bins.end(), find_fps),
                               result_user_bins.end());

        // insert false negatives
        // currently we don't have false negatives
        if (!fns_truths_per_query[id].empty())
            std::cerr << "warning: there is a FN I want to insert.\n";
        // auto & list = fns_truths_per_query[id];
        // result_user_bins.insert(result_user_bins.end(), list.begin(), list.end());
        // std::sort(result_user_bins.begin(), result_user_bins.end()); // sort again.

        if (result_user_bins.empty())
        {
            raptor_result_out << id << '\t' << '\n';
            continue;
        }

        // Write new normalised raptor line:
        auto it = result_user_bins.begin();
        raptor_result_out << id << '\t' << *(it++);
        while (it != result_user_bins.end())
            raptor_result_out << ',' << *(it++);
        raptor_result_out << '\n';
    }
    while (std::getline(raptor_result_in, line_buffer));

    std::cerr << "Done" << std::endl;
}

void init_parser(seqan3::argument_parser & parser, config & cfg)
{
    parser.add_option(cfg.raptor_result_file,
                      '\0',
                      "raptor_results",
                      "The raptor result file, e.g., \"raptor.results\".",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{});

    parser.add_option(cfg.FN_file,
                      '\0',
                      "fns",
                      "The true false negatives to incooporate.\".",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{});

    parser.add_option(cfg.FP_file,
                      '\0',
                      "fps",
                      "The true false negatives to incooporate.\".",
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
    seqan3::argument_parser parser{"normalise_mantis_output", argc, argv, seqan3::update_notifications::off};
    parser.info.author = "Svenja Mehringer, Enrico Seiler";
    parser.info.email = "enrico.seiler@fu-berlin.de";
    parser.info.short_description = "Corrects the raptor file with validated FPs and FNs.";
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

    correct_truth_file(cfg);
}

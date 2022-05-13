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

#include "check_output_file.hpp"
#include "parse_user_bin_ids.hpp"

struct config
{
    std::filesystem::path truth_user_bin_ids_file{};
    std::filesystem::path raptor_result_file{};
    std::filesystem::path output_file{};
};

robin_hood::unordered_map<std::string, uint64_t> create_ub_to_ub_mapping_from_header(std::ifstream & raptor_result_in,
                                                                                     std::string & line_buffer,
                                                                                     config const & cfg)
{
    // map[reference_name] = number
    std::cerr << "Reading " << cfg.truth_user_bin_ids_file << " ... " << std::flush;
    robin_hood::unordered_map<std::string, uint64_t> const truth_ub_name_to_id{parse_user_bin_ids(cfg.truth_user_bin_ids_file)};
    std::cerr << "Done" << std::endl;

    robin_hood::unordered_map<std::string, uint64_t> ub_to_ub;

    std::cerr << "Create ub_to_ub mapping ... "  << std::flush;
    while (std::getline(raptor_result_in, line_buffer) && line_buffer.starts_with('#') && line_buffer[1] != ('Q'))
    {
        auto tab_it{line_buffer.begin() + line_buffer.find('\t')};
        std::string_view const idx{line_buffer.begin() + 1/* skip '#' */,
                                   line_buffer.begin() + line_buffer.find('\t')};
        std::string_view const name_key{line_buffer.begin() + line_buffer.find_last_of('/') + 1,
                                        line_buffer.begin() + line_buffer.find(".fna.gz")};

        // std::cerr << "SEARCH FOR" << name_key << std::endl;

        ub_to_ub.emplace(idx, truth_ub_name_to_id.at(std::string{name_key}));
    }
    std::cerr << "Done" << std::endl;

    return ub_to_ub;
}

void normalise_output(config const & cfg)
{
    // Process raptor results
    std::ifstream raptor_result_in{cfg.raptor_result_file};
    std::ofstream raptor_result_out{cfg.output_file};

    // Buffers for file I/O
    std::vector<uint64_t> result_user_bins{};
    std::string line_buffer{};

    // ## Raptor results ##
    // The header stores the user bin ID.
    // The ids of the file to normalize have to be adapted to the truth sets user bin ids
    auto const ub_to_ub = create_ub_to_ub_mapping_from_header(raptor_result_in, line_buffer, cfg);

    std::cerr << "Processing " << cfg.raptor_result_file << " ... " << std::flush;

    while (std::getline(raptor_result_in, line_buffer))
    {
        result_user_bins.clear();

        auto tab_it{line_buffer.begin() + line_buffer.find('\t')};
        std::string_view const id{line_buffer.begin(), tab_it};
        std::string_view const bins{++tab_it, line_buffer.end()};

        std::string_view::size_type current_pos = 0;
        std::string_view::size_type comma_pos{bins.find(',')};

        while (comma_pos != std::string_view::npos)
        {
            result_user_bins.push_back(ub_to_ub.at(std::string(&bins[current_pos], comma_pos - current_pos)));
            current_pos = comma_pos + 1;
            comma_pos = bins.find(',', current_pos);
        }
        result_user_bins.push_back(ub_to_ub.at(std::string(&bins[current_pos], bins.size() - current_pos)));
        std::sort(result_user_bins.begin(), result_user_bins.end()); // compare script afterwards requires sorted UBs

        // Write new normalised raptor line:
        auto it = result_user_bins.begin();
        raptor_result_out << id << '\t' << *(it++);
        while (it != result_user_bins.end())
            raptor_result_out << ',' << *(it++);
        raptor_result_out << '\n';
    }

    std::cerr << "Done" << std::endl;
}

void init_parser(seqan3::argument_parser & parser, config & cfg)
{
    parser.add_option(cfg.truth_user_bin_ids_file,
                      '\0',
                      "user_bin_ids",
                      "The file containing user bin ids from the 'truth'-raptor file, e.g., \"user_bin.ids\".",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{});
    parser.add_option(cfg.raptor_result_file,
                      '\0',
                      "raptor_results",
                      "The raptor result file, e.g., \"raptor.results\".",
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
    parser.info.short_description = "Unifies raptor results by replacing user bin ids from one raptor file with those "
                                    "of a 'truth'-raptor file.";
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

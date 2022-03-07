// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <cassert>
#include <charconv>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <ranges>
#include <string>
#include <vector>

#include <seqan3/argument_parser/argument_parser.hpp>

#include "check_output_file.hpp"
#include "parse_user_bin_ids.hpp"

struct config
{
    std::filesystem::path mantis_result_file{};
    std::filesystem::path raptor_result_file{};
    std::filesystem::path user_bin_ids_file{};
    std::filesystem::path output_directory{};
};

auto find_tab(std::string const & str)
{
    auto const pos = str.find('\t');
    if (pos == std::string::npos)
        throw std::runtime_error{"Line \"" + str + "\" does not contain a tab."};
    return str.begin() + pos;
}

auto extract_hit(auto const & outer_it)
{
    auto const inner_it = *outer_it;
    char const * const first_char = std::addressof(*(inner_it.begin()));
    auto const length = std::ranges::distance(inner_it);

    uint64_t result{};
    std::from_chars(first_char, first_char + length, result);
    return result;
}

void compare_results(config const & cfg)
{
    // map[reference_name] = number
    robin_hood::unordered_map<std::string, uint64_t> const ub_name_to_id{parse_user_bin_ids(cfg.user_bin_ids_file)};

    std::ifstream mantis_result{cfg.mantis_result_file};
    std::ifstream raptor_result{cfg.raptor_result_file};

    std::ofstream false_positives_file{cfg.output_directory / "raptor.fps"};
    std::ofstream false_negatives_file{cfg.output_directory / "raptor.fns"};
    std::ofstream missing_ground_truths_file{cfg.output_directory / "missing_ground_truth.warn"};
    std::ofstream missing_lines_file{cfg.output_directory / "missing_lines.warn"};
    std::ofstream stats_file{cfg.output_directory / "stats.tsv"};

    bool missing_lines{false};

    std::string mantis_line{};
    std::string raptor_line{};
    uint64_t false_positives{};
    uint64_t false_negatives{};
    uint64_t line_no{};
    uint64_t mantis_hit_count{};
    uint64_t mantis_miss{};
    uint64_t raptor_hit_count{};
    uint64_t raptor_miss{};

    std::string query_name_buffer{};
    auto parse_original_bin = [&query_name_buffer, &ub_name_to_id] (std::string_view const & line)
    {
        // E.g., "GCF_000005825.2_ASM582v2_genomic106". 106 is the read number.
        // find() returns an iterator to the "g" of "genomic". `+7` moves the iterator to "1", which is the end
        // of the bin name.
        query_name_buffer.assign(line.begin(), line.begin() + line.find("genomic") + 7);
        return ub_name_to_id.at(query_name_buffer);
    };

    auto parse_query_name = [&mantis_line, &raptor_line] (auto const & mantis_tab_it, auto const & raptor_tab_it)
    {
        std::string_view const mantis_query_name{mantis_line.begin(), mantis_tab_it};
        std::string_view const raptor_query_name{raptor_line.begin(), raptor_tab_it};
        if (mantis_query_name != raptor_query_name)
            throw std::runtime_error{"Query names do not match, something went wrong"};
        return mantis_query_name;
    };

    while (std::getline(mantis_result, mantis_line) && std::getline(raptor_result, raptor_line))
    {
        auto const mantis_tab_it{find_tab(mantis_line)};
        auto const raptor_tab_it{find_tab(raptor_line)};
        std::string_view const query_name{parse_query_name(mantis_tab_it, raptor_tab_it)};

        uint64_t const original_bin{parse_original_bin(query_name)};
        bool mantis_found_correct_bin{false};
        bool raptor_found_correct_bin{false};

        std::ranges::split_view const mantis_fields_view{std::string_view{mantis_tab_it + 1, mantis_line.end()}, ','};
        std::ranges::split_view const raptor_fields_view{std::string_view{raptor_tab_it + 1, raptor_line.end()}, ','};
        auto mantis_it{mantis_fields_view.begin()};
        auto raptor_it{raptor_fields_view.begin()};

        while (mantis_it != mantis_fields_view.end() && raptor_it != raptor_fields_view.end())
        {
            uint64_t const mantis_hit_bin{extract_hit(mantis_it)};
            mantis_found_correct_bin = mantis_found_correct_bin || mantis_hit_bin == original_bin;

            uint64_t const raptor_hit_bin{extract_hit(raptor_it)};
            raptor_found_correct_bin = raptor_found_correct_bin || raptor_hit_bin == original_bin;

            if (mantis_hit_bin != raptor_hit_bin)
            {
                if (mantis_hit_bin < raptor_hit_bin)
                {
                    if (raptor_hit_bin != original_bin)
                    {
                        false_negatives_file << query_name << ":" << mantis_hit_bin << '\n';
                        ++false_negatives;
                    }
                    ++mantis_hit_count;
                    ++mantis_it;
                }
                else
                {
                    if (raptor_hit_bin != original_bin)
                    {
                        false_positives_file << query_name << ":" << raptor_hit_bin << '\n';
                        ++false_positives;
                    }
                    ++raptor_hit_count;
                    ++raptor_it;
                }
            }
            else
            {
                 ++mantis_hit_count;
                 ++raptor_hit_count;
                 ++mantis_it;
                 ++raptor_it;
            }
        }

        // process the rest of mantis
        while (mantis_it != mantis_fields_view.end())
        {
            uint64_t const mantis_hit_bin{extract_hit(mantis_it)};
            mantis_found_correct_bin = mantis_found_correct_bin || mantis_hit_bin == original_bin;
            false_negatives_file << query_name << ":" << mantis_hit_bin << '\n';
            ++false_negatives;
            ++mantis_hit_count;
            ++mantis_it;
        }

        // process the rest of raptor
        while (raptor_it != raptor_fields_view.end())
        {
            uint64_t const raptor_hit_bin{extract_hit(raptor_it)};
            raptor_found_correct_bin = raptor_found_correct_bin || raptor_hit_bin == original_bin;
            false_positives_file << query_name << ":" << raptor_hit_bin << '\n';
            ++false_positives;
            ++raptor_hit_count;
            ++raptor_it;
        }

        if (!mantis_found_correct_bin)
        {
            ++mantis_miss;
            missing_ground_truths_file << "Line " << line_no << ": "
                                       << "Could not find query " << query_name << " "
                                       << "(" << query_name_buffer << ":" << original_bin << ") "
                                       << "in its respective genome in mantis.\n";
        }
        if (!raptor_found_correct_bin)
        {
            ++raptor_miss;
            missing_ground_truths_file << "Line " << line_no << ": "
                                       << "Could not find query " << query_name << " "
                                       << "(" << query_name_buffer << ":" << original_bin << ") "
                                       << "in its respective genome in raptor.\n";
        }

        ++line_no;
    }

    while (std::getline(mantis_result, mantis_line))
    {
        missing_lines = true;
        missing_lines_file << "Missing line of mantis in comparison: " << mantis_line << '\n';
    }
    while (std::getline(raptor_result, raptor_line))
    {
        missing_lines = true;
        missing_lines_file << "Missing line of raptor in comparison: " << raptor_line << '\n';
    }

    stats_file << "Mantis total:\t" << mantis_hit_count << '\n';
    stats_file << "Mantis miss:\t" << mantis_miss << '\n';
    stats_file << "Raptor total:\t" << raptor_hit_count << '\n';
    stats_file << "Raptor miss:\t" << raptor_miss << '\n';
    stats_file << "Raptor FP:\t" << false_positives << '\n';
    stats_file << "Raptor FN:\t" << false_negatives << '\n';

    if (missing_lines)
        std::cout << "[WARNING] Somes lines were missing. See "
                  << cfg.output_directory.c_str()
                  << "missing_lines.warn\n";

    if (mantis_miss || raptor_miss)
        std::cout << "[Info] Missing ground truths are listed in "
                  << cfg.output_directory.c_str()
                  << "missing_ground_truth.warn\n";

    std::cout << "[Info] False positives: "
                << cfg.output_directory.c_str()
                << "raptor.fps\n";

    std::cout << "[Info] False negatives: "
                << cfg.output_directory.c_str()
                << "raptor.fns\n";

    std::cout << "[Info] Statistics: "
                << cfg.output_directory.c_str()
                << "stats.tsv\n";

    std::cout << "[Info] Content of stats.tsv:\n"
              << "       Mantis total:\t" << mantis_hit_count << '\n'
              << "       Mantis miss:\t" << mantis_miss << '\n'
              << "       Raptor total:\t" << raptor_hit_count << '\n'
              << "       Raptor miss:\t" << raptor_miss << '\n'
              << "       Raptor FP:\t" << false_positives << '\n'
              << "       Raptor FN:\t" << false_negatives << '\n';
}

void init_parser(seqan3::argument_parser & parser, config & cfg)
{
    parser.add_option(cfg.mantis_result_file,
                      '\0',
                      "mantis_results",
                      "The mantis result file produced by normalise_mantis_output.",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{});
    parser.add_option(cfg.raptor_result_file,
                      '\0',
                      "raptor_results",
                      "The raptor result file, e.g., \"raptor.results\".",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{});
    parser.add_option(cfg.user_bin_ids_file,
                      '\0',
                      "user_bin_ids",
                      "The file containing user bin ids, e.g., \"user_bin.ids\".",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{});
    parser.add_option(cfg.output_directory,
                      '\0',
                      "output_directory",
                      "Provide a path to the output.",
                      seqan3::option_spec::required);
}

int main(int argc, char ** argv)
{
    seqan3::argument_parser parser{"compare_mantis_raptor_output", argc, argv, seqan3::update_notifications::off};
    parser.info.author = "Svenja Mehringer, Enrico Seiler";
    parser.info.email = "enrico.seiler@fu-berlin.de";
    parser.info.short_description = "Compares mantis and raptor results.";
    parser.info.version = "0.0.1";

    config cfg{};
    init_parser(parser, cfg);

    try
    {
        parser.parse();
        cfg.output_directory = std::filesystem::absolute(cfg.output_directory);
        check_output_file(cfg.output_directory / "stats.txt");
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        std::cerr << "[Error] " << ext.what() << '\n';
        std::exit(-1);
    }

    compare_results(cfg);
}

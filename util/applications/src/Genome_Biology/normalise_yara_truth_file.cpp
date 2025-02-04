// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <cassert>
#include <charconv>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <sharg/parser.hpp>

#include <hibf/contrib/robin_hood.hpp>

inline void check_output_file(std::filesystem::path const & output_file)
{
    std::filesystem::path const output_directory = output_file.parent_path();
    std::error_code ec{};
    std::filesystem::create_directories(output_directory, ec);

    if (!output_directory.empty() && ec)
        sharg::parser_error{
            sharg::detail::to_string("Failed to create directory\"", output_directory.c_str(), "\": ", ec.message())};
}

struct options
{
    std::filesystem::path yara_result_file{};
    std::filesystem::path output_file{};
    std::filesystem::path query_names_file{};
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
    std::cerr << "Done" << '\n';
    return query_names;
}

void normalise_output(options const & cfg)
{
    std::vector<std::string> const query_names{parse_query_names(cfg.query_names_file)};
    std::cerr << "Read " << query_names.size() << "query names" << '\n';

    // Process yara results
    std::ifstream yara_result_in{cfg.yara_result_file};
    std::ofstream yara_result_out{cfg.output_file};
    std::vector<uint64_t> results;
    std::string last_seen_query_name{};

    // Buffers for file I/O
    std::string result_buffer{};
    std::string line_buffer{};

    auto parse_query_name_and_user_bin = [](std::string const & line)
    {
        uint64_t idx{};
        std::string_view const qname{line.begin(), line.begin() + line.find(':')};
        std::string_view const idx_str{line.begin() + qname.size() + 1, line.end()};
        std::from_chars(idx_str.data(), idx_str.data() + idx_str.size(), idx);

        return std::make_pair(qname, idx);
    };

    auto process_results = [&results, &result_buffer, &yara_result_out]()
    {
        if (!results.empty())
        {
            std::ranges::sort(results);
            for (size_t const ub : results)
                result_buffer += std::to_string(ub) + ',';
            result_buffer.back() = '\n';

            yara_result_out << result_buffer;
            result_buffer.clear();
            results.clear();
        }
        else
        {
            yara_result_out << '\n';
        }
    };

    auto qname_it = query_names.begin();

    auto check_qname_against_reference = [&qname_it, &query_names, &last_seen_query_name, &yara_result_out]()
    {
        if (qname_it == query_names.end())
            throw std::runtime_error{"query_names consumed although processing has not ended. last_seen: "
                                     + last_seen_query_name};

        while (qname_it != query_names.end() && *qname_it != last_seen_query_name)
        {
            std::cerr << "Note: " << *qname_it << " not found in validation file." << '\n';
            yara_result_out << *qname_it << '\t' << '\n';
            ++qname_it;
        }

        if (qname_it != query_names.end())
            ++qname_it; // function is called when last_seen_query_name is also updated
    };

    std::cerr << "Processing " << cfg.yara_result_file << " ... " << std::flush;

    // First line.
    if (std::getline(yara_result_in, line_buffer))
    {
        auto && [qname, user_bin_idx] = parse_query_name_and_user_bin(line_buffer);
        last_seen_query_name = qname;
        results.push_back(user_bin_idx);
    }

    while (std::getline(yara_result_in, line_buffer))
    {
        auto && [qname, user_bin_idx] = parse_query_name_and_user_bin(line_buffer);

        if (qname != last_seen_query_name) // new query
        {
            check_qname_against_reference();

            yara_result_out << last_seen_query_name << '\t';

            process_results();

            results.clear();
            last_seen_query_name = qname;
            results.push_back(user_bin_idx);
        }
        else
        {
            results.push_back(user_bin_idx);
        }
    }

    // Write last results.
    yara_result_out << last_seen_query_name << '\t';
    process_results();
    check_qname_against_reference();

    if (qname_it != query_names.end())
        throw std::runtime_error{"query_names not fully consumed although processing has ended. last qname: "
                                 + (*qname_it)};

    std::cerr << "Done" << '\n';
}

void init_parser(sharg::parser & parser, options & cfg)
{
    parser.add_option(cfg.yara_result_file,
                      sharg::config{.short_id = '\0',
                                    .long_id = "yara_results",
                                    .description = "The yara result file, e.g., \"yara.results\".",
                                    .required = true});
    parser.add_option(cfg.query_names_file,
                      sharg::config{.short_id = '\0',
                                    .long_id = "query_names",
                                    .description = "The file containing query names, e.g., \"query.names\".",
                                    .required = true});
    parser.add_option(cfg.output_file,
                      sharg::config{.short_id = '\0',
                                    .long_id = "output_file",
                                    .description = "Provide a path to the output.",
                                    .required = true});
}

int main(int argc, char ** argv)
{
    sharg::parser parser{"normalise_yara_output", argc, argv, sharg::update_notifications::off};
    parser.info.author = "Svenja Mehringer, Enrico Seiler";
    parser.info.email = "enrico.seiler@fu-berlin.de";
    parser.info.short_description = "Converts yara results into raptor-like results.";
    parser.info.version = "0.0.1";

    options cfg{};
    init_parser(parser, cfg);

    try
    {
        parser.parse();
        check_output_file(cfg.output_file);
    }
    catch (sharg::parser_error const & ext)
    {
        std::cerr << "[Error] " << ext.what() << '\n';
        std::exit(-1);
    }

    normalise_output(cfg);
}

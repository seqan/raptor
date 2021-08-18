// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <raptor/argument_parsing/build.hpp>
#include <raptor/build/build.hpp>

namespace raptor
{

void init_build_parser(seqan3::argument_parser & parser, build_arguments & arguments)
{
    init_shared_meta(parser);
    init_shared_options(parser, arguments);
    parser.add_positional_option(arguments.bin_file,
                                 arguments.is_socks ? "File containing color and file names." :
                                                      "File containing one file per line per bin.",
                                 seqan3::input_file_validator{});
    parser.add_option(arguments.parts,
                      '\0',
                      "parts",
                      "Splits the index in this many parts.",
                      arguments.is_socks ? seqan3::option_spec::hidden : seqan3::option_spec::standard,
                      power_of_two_validator{});
    parser.add_option(arguments.window_size,
                      '\0',
                      "window",
                      "Choose the window size.",
                      arguments.is_socks ? seqan3::option_spec::hidden : seqan3::option_spec::standard,
                      positive_integer_validator{});
    parser.add_option(arguments.kmer_size,
                      '\0',
                      "kmer",
                      "Choose the kmer size.",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{1, 32});
    parser.add_option(arguments.out_path,
                      '\0',
                      "output",
                      arguments.is_socks ? "Provide an output filepath." :
                                           "Provide an output filepath or an output directory if --compute-minimiser is used.",
                      seqan3::option_spec::required);
    parser.add_option(arguments.size,
                      '\0',
                      "size",
                      "Choose the size of the resulting index.",
                      seqan3::option_spec::required,
                      size_validator{"\\d+\\s{0,1}[k,m,g,t,K,M,G,T]"});
    parser.add_option(arguments.hash,
                      '\0',
                      "hash",
                      "Choose the number of hashes.",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{1, 5});
    parser.add_flag(arguments.compressed,
                    '\0',
                    "compressed",
                    "Build a compressed index.");
    parser.add_flag(arguments.compute_minimiser,
                    '\0',
                    "compute-minimiser",
                    "Computes minimisers using cutoffs from Mantis (Pandey et al.). Does not create the index.",
                    arguments.is_socks ? seqan3::option_spec::hidden : seqan3::option_spec::standard);
}

void run_build(seqan3::argument_parser & parser, bool const is_socks)
{
    build_arguments arguments{};
    arguments.is_socks = is_socks;
    init_build_parser(parser, arguments);
    try_parsing(parser);

    // ==========================================
    // Process bin_path
    // ==========================================
    if (!arguments.is_socks) // Either only one bin or a file containing bin paths
    {
        std::ifstream istrm{arguments.bin_file};
        std::string line;
        auto sequence_file_validator{bin_validator{}.sequence_file_validator};

        while (std::getline(istrm, line))
        {
            if (!line.empty())
            {
                sequence_file_validator(line);
                arguments.bin_path.emplace_back(std::vector<std::string>{line});
            }
        }
    }
    else
    {
        std::ifstream istrm{arguments.bin_file};
        std::string line;
        std::string color_name;
        std::string file_name;
        std::vector<std::string> tmp;
        auto sequence_file_validator{bin_validator{}.sequence_file_validator};

        while (std::getline(istrm, line))
        {
            if (!line.empty())
            {
                tmp.clear();
                std::stringstream sstream{line};
                sstream >> color_name;
                while (std::getline(sstream, file_name, ' '))
                {
                    if (!file_name.empty())
                    {
                        sequence_file_validator(file_name);
                        tmp.emplace_back(file_name);
                    }
                }
                arguments.bin_path.emplace_back(tmp);
            }
        }
    }

    // ==========================================
    // Various checks.
    // ==========================================

    arguments.bins = arguments.bin_path.size();

    if (parser.is_option_set("window"))
    {
        if (arguments.kmer_size > arguments.window_size)
            throw seqan3::argument_parser_error{"The k-mer size cannot be bigger than the window size."};
    }
    else
        arguments.window_size = arguments.kmer_size;

    if (parser.is_option_set("compute-minimiser"))
    {
        try
        {
            seqan3::output_directory_validator{}(arguments.out_path);
        }
        catch (seqan3::argument_parser_error const & ext)
        {
            std::cerr << "[Error] " << ext.what() << '\n';
            std::exit(-1);
        }
    }
    else
    {
        try
        {
            seqan3::output_file_validator{}(arguments.out_path);
        }
        catch (seqan3::argument_parser_error const & ext)
        {
            std::cerr << "[Error] " << ext.what() << '\n';
            std::exit(-1);
        }
    }

    // ==========================================
    // Process --size.
    // ==========================================
    arguments.size.erase(std::remove(arguments.size.begin(), arguments.size.end(), ' '), arguments.size.end());

    size_t multiplier{};

    switch (std::tolower(arguments.size.back()))
    {
        case 't':
            multiplier = 8ull * 1024ull * 1024ull * 1024ull * 1024ull;
            break;
        case 'g':
            multiplier = 8ull * 1024ull * 1024ull * 1024ull;
            break;
        case 'm':
            multiplier = 8ull * 1024ull * 1024ull;
            break;
        case 'k':
            multiplier = 8ull * 1024ull;
            break;
        default:
            throw seqan3::argument_parser_error{"Use {k, m, g, t} to pass size. E.g., --size 8g."};
    }

    size_t size{};
    std::from_chars(arguments.size.data(), arguments.size.data() + arguments.size.size() - 1, size);
    size *= multiplier;
    arguments.bits = size / (((arguments.bins + 63) >> 6) << 6);

    // ==========================================
    // Read w and k from minimiser header file
    // ==========================================
    if (std::filesystem::path header_file_path = arguments.bin_path[0][0]; header_file_path.extension() == ".minimiser")
    {
        header_file_path.replace_extension("header");
        std::ifstream file_stream{header_file_path};
        uint64_t kmer_size{};
        file_stream >> kmer_size >> arguments.window_size;
        arguments.kmer_size = kmer_size;
    }

    // ==========================================
    // Dispatch
    // ==========================================
    raptor_build(arguments);
};

} // namespace raptor

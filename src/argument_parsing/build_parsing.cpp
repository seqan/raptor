// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <raptor/argument_parsing/build_parsing.hpp>
#include <raptor/argument_parsing/init_shared_meta.hpp>
#include <raptor/argument_parsing/parse_bin_path.hpp>
#include <raptor/argument_parsing/validators.hpp>
#include <raptor/build/raptor_build.hpp>

namespace raptor
{

void init_build_parser(seqan3::argument_parser & parser, build_arguments & arguments)
{
    init_shared_meta(parser);
    parser.info.examples = {"raptor build --kmer 19 --window 23 --size 8m --output raptor.index all_bin_paths.txt",
                            "raptor build --kmer 19 --window 23 --compute-minimiser --output precomputed_minimisers all_bin_paths.txt",
                            "raptor build --size 8m --output minimiser_raptor.index all_minimiser_paths.txt"};
    parser.add_positional_option(arguments.bin_file,
                                 (arguments.is_socks ? "File containing color and file names. " :
                                                       "File containing file names. ") + bin_validator{}.get_help_page_message(),
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
                      "The window size. \\fI\\fBIf not set, defaults to the k-mer size.\\fP",
                      arguments.is_socks ? seqan3::option_spec::hidden : seqan3::option_spec::standard,
                      positive_integer_validator{});
    parser.add_option(arguments.kmer_size,
                      '\0',
                      "kmer",
                      "The k-mer size.", // Mutually exclusive with --shape.
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{1, 32});
    parser.add_option(arguments.shape_string,
                      '\0',
                      "shape",
                      "The shape to use for k-mers. Mutually exclusive with --kmer.",
                      seqan3::option_spec::advanced, // Add help in kmer_size
                      seqan3::regex_validator{"[01]+"});
    parser.add_option(arguments.out_path,
                      '\0',
                      "output",
                      arguments.is_socks ? "Provide an output filepath." :
                                           "Provide an output filepath or an output directory if --compute-minimiser is used.",
                      seqan3::option_spec::required);
    parser.add_option(arguments.size,
                      '\0',
                      "size",
                      "The size in bytes of the resulting index.",
                      seqan3::option_spec::standard,
                      size_validator{"\\d+\\s{0,1}[k,m,g,t,K,M,G,T]"});
    parser.add_option(arguments.hash,
                      '\0',
                      "hash",
                      "The number of hash functions to use.",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{1, 5});
    parser.add_option(arguments.threads,
                      '\0',
                      "threads",
                      "The numer of threads to use.",
                      seqan3::option_spec::standard,
                      positive_integer_validator{});
    parser.add_option(arguments.fpr,
                      '\0',
                      "fpr",
                      "False positive rate of the HIBF.",
                      seqan3::option_spec::advanced,
                      seqan3::arithmetic_range_validator{0.0, 1.0});
    parser.add_flag(arguments.compressed,
                    '\0',
                    "compressed",
                    "Build a compressed index.");
    parser.add_flag(arguments.compute_minimiser,
                    '\0',
                    "compute-minimiser",
                    "Computes minimisers using cutoffs from Mantis (Pandey et al.). Does not create the index.",
                    arguments.is_socks ? seqan3::option_spec::hidden : seqan3::option_spec::standard);
    parser.add_flag(arguments.compute_minimiser,
                    '\0',
                    "compute-minimizer",
                    "Hidden flag, alias of --compute-minimiser.",
                    seqan3::option_spec::hidden);
    parser.add_flag(arguments.disable_cutoffs,
                    '\0',
                    "disable-cutoffs",
                    "Do not apply cutoffs when using --compute-minimiser.",
                    arguments.is_socks ? seqan3::option_spec::hidden : seqan3::option_spec::standard);
    parser.add_flag(arguments.is_hibf,
                    '\0',
                    "hibf",
                    "Index is an HIBF.",
                    seqan3::option_spec::advanced);
}

void build_parsing(seqan3::argument_parser & parser, bool const is_socks)
{
    build_arguments arguments{};
    arguments.is_socks = is_socks;
    init_build_parser(parser, arguments);
    parser.parse();

    // ==========================================
    // Various checks.
    // ==========================================
    if (parser.is_option_set("shape"))
    {
        if (parser.is_option_set("kmer"))
            throw seqan3::argument_parser_error{"You cannot set both shape and k-mer arguments."};

        uint64_t tmp{};

        std::from_chars(arguments.shape_string.data(),
                        arguments.shape_string.data() + arguments.shape_string.size(),
                        tmp,
                        2);

        arguments.shape = seqan3::shape{seqan3::bin_literal{tmp}};
    }
    else
    {
        arguments.shape = seqan3::shape{seqan3::ungapped{arguments.kmer_size}};
    }

    if (parser.is_option_set("window"))
    {
        if (arguments.shape.size() > arguments.window_size)
            throw seqan3::argument_parser_error{"The k-mer size cannot be bigger than the window size."};
    }
    else
    {
        arguments.window_size = arguments.shape.size();
    }

    bool const is_compute_minimiser_set{parser.is_option_set("compute-minimiser") ||
                                        parser.is_option_set("compute-minimizer")};

    arguments.compute_minimiser = is_compute_minimiser_set;

    std::filesystem::path output_directory = is_compute_minimiser_set ? arguments.out_path :
                                                                        arguments.out_path.parent_path();
    std::error_code ec{};
    std::filesystem::create_directories(output_directory, ec);

// LCOV_EXCL_START
    if (!output_directory.empty() && ec)
        throw seqan3::argument_parser_error{seqan3::detail::to_string("Failed to create directory\"",
                                                                      output_directory.c_str(),
                                                                      "\": ",
                                                                      ec.message())};
// LCOV_EXCL_STOP

    if (!is_compute_minimiser_set)
    {
        seqan3::output_file_validator{}(arguments.out_path);

        if (!parser.is_option_set("size") && !parser.is_option_set("hibf"))
        {
            throw seqan3::argument_parser_error{"Option --size is required but not set."};
        }
    }

    // ==========================================
    // Process bin_path
    // ==========================================
    parse_bin_path(arguments);
    arguments.bins = arguments.bin_path.size();

    // ==========================================
    // Process --size.
    // ==========================================
    arguments.size.erase(std::remove(arguments.size.begin(), arguments.size.end(), ' '), arguments.size.end());

    size_t multiplier{};

    switch (std::tolower(arguments.size.back()))
    {
// LCOV_EXCL_START
        case 't':
            multiplier = 8ull * 1024ull * 1024ull * 1024ull * 1024ull;
            break;
        case 'g':
            multiplier = 8ull * 1024ull * 1024ull * 1024ull;
            break;
        case 'm':
            multiplier = 8ull * 1024ull * 1024ull;
            break;
// LCOV_EXCL_STOP
        case 'k':
            multiplier = 8ull * 1024ull;
            break;
// LCOV_EXCL_START
        default:
            throw seqan3::argument_parser_error{"Use {k, m, g, t} to pass size. E.g., --size 8g."};
// LCOV_EXCL_STOP
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
        std::string shape_string{};
        file_stream >> shape_string >> arguments.window_size;

        uint64_t tmp{};
        std::from_chars(shape_string.data(),
                        shape_string.data() + shape_string.size(),
                        tmp,
                        2);

        arguments.shape = seqan3::shape{seqan3::bin_literal{tmp}};
    }

    // ==========================================
    // Dispatch
    // ==========================================
    raptor_build(arguments);
};

} // namespace raptor

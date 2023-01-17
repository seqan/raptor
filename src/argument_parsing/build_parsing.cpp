// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <raptor/argument_parsing/build_parsing.hpp>
#include <raptor/argument_parsing/init_shared_meta.hpp>
#include <raptor/argument_parsing/parse_bin_path.hpp>
#include <raptor/argument_parsing/validators.hpp>
#include <raptor/build/raptor_build.hpp>

namespace raptor
{

void init_build_parser(sharg::parser & parser, build_arguments & arguments)
{
    init_shared_meta(parser);
    parser.info.description.emplace_back("Constructs a Raptor index of large collections of nucleotide sequences.");
    parser.info.examples = {
        "raptor build --kmer 19 --window 23 --size 8m --output raptor.index all_bin_paths.txt",
        "raptor build --kmer 19 --window 23 --compute-minimiser --output precomputed_minimisers all_bin_paths.txt",
        "raptor build --size 8m --output minimiser_raptor.index all_minimiser_paths.txt"};

    parser.add_positional_option(
        arguments.bin_file,
        sharg::config{.description = "File containing file names. " + bin_validator{}.get_help_page_message(),
                      .validator = sharg::input_file_validator{}});

    parser.add_option(arguments.parts,
                      sharg::config{.short_id = '\0',
                                    .long_id = "parts",
                                    .description = "Splits the index in this many parts.",
                                    .validator = power_of_two_validator{}});
    parser.add_option(arguments.kmer_size,
                      sharg::config{.short_id = '\0',
                                    .long_id = "kmer",
                                    .description = "The k-mer size.",
                                    .validator = sharg::arithmetic_range_validator{1, 32}});
    parser.add_option(arguments.window_size,
                      sharg::config{.short_id = '\0',
                                    .long_id = "window",
                                    .description = "The window size.",
                                    .default_message = "k-mer size",
                                    .validator = positive_integer_validator{}});
    parser.add_option(arguments.shape_string,
                      sharg::config{.short_id = '\0',
                                    .long_id = "shape",
                                    .description = "The shape to use for k-mers. Mutually exclusive with --kmer.",
                                    .advanced = true,
                                    .validator = sharg::regex_validator{"[01]+"}});
    parser.add_option(
        arguments.out_path,
        sharg::config{.short_id = '\0',
                      .long_id = "output",
                      .description =
                          "Provide an output filepath or an output directory if --compute-minimiser is used.",
                      .required = true});
    parser.add_option(arguments.size,
                      sharg::config{.short_id = '\0',
                                    .long_id = "size",
                                    .description = "The size in bytes of the resulting index.",
                                    .validator = size_validator{"\\d+\\s{0,1}[k,m,g,t,K,M,G,T]"}});
    parser.add_option(arguments.hash,
                      sharg::config{.short_id = '\0',
                                    .long_id = "hash",
                                    .description = "The number of hash functions to use.",
                                    .validator = sharg::arithmetic_range_validator{1, 5}});
    parser.add_option(arguments.threads,
                      sharg::config{.short_id = '\0',
                                    .long_id = "threads",
                                    .description = "The number of threads to use.",
                                    .validator = positive_integer_validator{}});
    parser.add_option(arguments.fpr,
                      sharg::config{.short_id = '\0',
                                    .long_id = "fpr",
                                    .description = "False positive rate of the HIBF.",
                                    .advanced = true,
                                    .validator = sharg::arithmetic_range_validator{0.0, 1.0}});

    parser.add_flag(
        arguments.compressed,
        sharg::config{.short_id = '\0', .long_id = "compressed", .description = "Build a compressed index."});
    parser.add_flag(
        arguments.compute_minimiser,
        sharg::config{.short_id = '\0',
                      .long_id = "compute-minimiser",
                      .description =
                          "Computes minimisers using cutoffs from Mantis (Pandey et al.). Does not create the index."});
    parser.add_flag(arguments.compute_minimiser,
                    sharg::config{.short_id = '\0',
                                  .long_id = "compute-minimizer",
                                  .description = "Hidden flag, alias of --compute-minimiser.",
                                  .hidden = true});
    parser.add_flag(arguments.disable_cutoffs,
                    sharg::config{.short_id = '\0',
                                  .long_id = "disable-cutoffs",
                                  .description = "Do not apply cutoffs when using --compute-minimiser."});
    parser.add_flag(
        arguments.is_hibf,
        sharg::config{.short_id = '\0', .long_id = "hibf", .description = "Index is an HIBF.", .advanced = true});
}

void build_parsing(sharg::parser & parser)
{
    build_arguments arguments{};
    init_build_parser(parser, arguments);
    parser.parse();

    // ==========================================
    // Various checks.
    // ==========================================
    if (parser.is_option_set("shape"))
    {
        if (parser.is_option_set("kmer"))
            throw sharg::parser_error{"You cannot set both shape and k-mer arguments."};

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

    if (!parser.is_option_set("window"))
        arguments.window_size = arguments.shape.size();
    else if (arguments.shape.size() > arguments.window_size)
        throw sharg::parser_error{"The k-mer size cannot be bigger than the window size."};

    bool const is_compute_minimiser_set{parser.is_option_set("compute-minimiser")
                                        || parser.is_option_set("compute-minimizer")};

    arguments.compute_minimiser = is_compute_minimiser_set;

    std::filesystem::path output_directory =
        is_compute_minimiser_set ? arguments.out_path : arguments.out_path.parent_path();
    std::error_code ec{};
    std::filesystem::create_directories(output_directory, ec);

    // GCOVR_EXCL_START
    if (!output_directory.empty() && ec)
        throw sharg::parser_error{
            seqan3::detail::to_string("Failed to create directory\"", output_directory.c_str(), "\": ", ec.message())};
    // GCOVR_EXCL_STOP

    if (!is_compute_minimiser_set)
    {
        sharg::output_file_validator{}(arguments.out_path);

        if (!parser.is_option_set("size") && !parser.is_option_set("hibf"))
        {
            throw sharg::parser_error{"Option --size is required but not set."};
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
    if (!parser.is_option_set("hibf"))
    {
        arguments.size.erase(std::remove(arguments.size.begin(), arguments.size.end(), ' '), arguments.size.end());
        size_t multiplier{};

        switch (std::tolower(arguments.size.back()))
        {
            // GCOVR_EXCL_START
        case 't':
            multiplier = 8ull * 1024ull * 1024ull * 1024ull * 1024ull;
            break;
        case 'g':
            multiplier = 8ull * 1024ull * 1024ull * 1024ull;
            break;
        case 'm':
            multiplier = 8ull * 1024ull * 1024ull;
            break;
            // GCOVR_EXCL_STOP
        case 'k':
            multiplier = 8ull * 1024ull;
            break;
            // GCOVR_EXCL_START
        default:
            throw sharg::parser_error{"Use {k, m, g, t} to pass size. E.g., --size 8g."};
            // GCOVR_EXCL_STOP
        }

        size_t size{};
        std::from_chars(arguments.size.data(), arguments.size.data() + arguments.size.size() - 1, size);
        size *= multiplier;
        arguments.bits = size / (((arguments.bins + 63) >> 6) << 6);
    }

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
        std::from_chars(shape_string.data(), shape_string.data() + shape_string.size(), tmp, 2);

        arguments.shape = seqan3::shape{seqan3::bin_literal{tmp}};
        arguments.is_minimiser = true;
    }

    // ==========================================
    // Dispatch
    // ==========================================
    raptor_build(arguments);
}

} // namespace raptor

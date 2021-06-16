// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/eseiler/minimizer_thresholds/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Random data generator.
 */

#include <seqan3/std/filesystem>
#include <random>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/output.hpp>

void run_program(std::filesystem::path const & out_directory,
                 size_t const reference_length,
                 size_t const number_of_queries,
                 size_t const query_length,
                 uint8_t const min_error,
                 uint8_t const max_error,
                 size_t const seed)
{
    std::mt19937_64 gen(seed);
    std::uniform_int_distribution<uint8_t> dis_dna(0, 3); // dna4
    std::uniform_int_distribution<size_t> dis_pos(0, reference_length - query_length - 1); // positions
    std::uniform_int_distribution<size_t> dis_err(0, query_length - 1); // errors


    // =========================================================================
    // Reference
    // =========================================================================
    std::vector<seqan3::dna4> reference_sequence;
    reference_sequence.reserve(reference_length);

    for (size_t i = 0; i < reference_length; ++i)
        reference_sequence.push_back(seqan3::dna4{}.assign_rank(dis_dna(gen)));

    {
        std::filesystem::path reference_filename = out_directory / "reference.fasta";
        seqan3::sequence_file_output reference_out{reference_filename};
        reference_out.emplace_back(reference_sequence, "reference_" + std::to_string(reference_length));
    }

    // =========================================================================
    // Queries with fixed number of errors
    // =========================================================================
    std::vector<seqan3::dna4> query_sequence;
    query_sequence.reserve(query_length);

    for (uint8_t errors = min_error; errors <= max_error; ++errors)
    {
        std::filesystem::path query_filename = out_directory / ("query_e" + std::to_string(errors) + ".fasta");
        {
            seqan3::sequence_file_output query_out{query_filename};
            for (size_t i = 0; i < number_of_queries; ++i)
            {
                size_t start = dis_pos(gen);
                query_sequence = reference_sequence |
                                 seqan3::views::slice(start, start + query_length) |
                                 seqan3::views::to<std::vector<seqan3::dna4>>;

                for (size_t j = 0; j < errors; ++j)
                {
                    size_t error_pos = dis_err(gen);
                    seqan3::dna4 new_base{};
                    new_base.assign_rank(dis_dna(gen));
                    while (new_base == query_sequence[error_pos])
                        new_base.assign_rank(dis_dna(gen));
                    query_sequence[error_pos] = new_base;
                }
                query_out.emplace_back(query_sequence, "query_" + std::to_string(i));
            }
        }
    }

    // =========================================================================
    // Random queries (to check for false positives)
    // =========================================================================
    {
        std::filesystem::path query_filename = out_directory / "query_random.fasta";
        seqan3::sequence_file_output query_out{query_filename};

        for (size_t i = 0; i < number_of_queries; ++i)
        {
            query_sequence.clear();

            for (size_t i = 0; i < query_length; ++i)
                query_sequence.push_back(seqan3::dna4{}.assign_rank(dis_dna(gen)));

            query_out.emplace_back(query_sequence, "query_" + std::to_string(i));
        }
    }
}

struct cmd_arguments
{
    std::filesystem::path output_file{};
    size_t reference_length{10'000'000};
    size_t number_of_queries{100'000};
    size_t query_length{150};
    uint8_t min_error{0};
    uint8_t max_error{3};
    size_t seed{0x7E82B6F280D4706B};
};

void initialize_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "Enrico Seiler";
    parser.info.short_description = "This programs creates random data for testing purposes.";
    parser.info.version = "1.0.0";

    parser.add_option(args.output_file, '\0', "out", "The output directory where the files will be located.",
                      seqan3::option_spec::required, seqan3::output_directory_validator{});
    parser.add_option(args.reference_length, '\0', "reference-size", "The length of the reference.",
                      seqan3::option_spec::standard, seqan3::arithmetic_range_validator{1, 1'000'000'000});
    parser.add_option(args.number_of_queries, '\0', "number-of-queries", "The number of queries.",
                      seqan3::option_spec::standard, seqan3::arithmetic_range_validator{1, 1'000'000'000});
    parser.add_option(args.query_length, '\0', "query_length", "The length of the queries.",
                      seqan3::option_spec::standard, seqan3::arithmetic_range_validator{1, 1'000'000});
    parser.add_option(args.min_error, '\0', "min_error", "The minimal number of errors.");
    parser.add_option(args.max_error, '\0', "max_error", "The maximal number of errors.");
    parser.add_option(args.seed, '\0', "seed", "The seed to use.", seqan3::option_spec::advanced);
}

int main(int argc, char ** argv)
{
    seqan3::argument_parser myparser{"Minimizers", argc, argv, seqan3::update_notifications::off};
    cmd_arguments args{};

    initialize_argument_parser(myparser, args);

    try
    {
        myparser.parse();
        if (args.min_error > args.max_error)
           throw seqan3::parser_invalid_argument{"The minimum number of errors cannot be greater than the maximum."};
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        std::cerr << "[Error] " << ext.what() << "\n";
        return -1;
    }

    run_program(args.output_file,
                args.reference_length,
                args.number_of_queries,
                args.query_length,
                args.min_error,
                args.max_error,
                args.seed);

    return 0;
}

// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <random>

#include <sharg/all.hpp>

#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>

struct my_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

struct cmd_arguments
{
    std::filesystem::path bin_file_path{};
    std::vector<std::filesystem::path> bin_path{};
    std::filesystem::path out_path{};
    // uint8_t min_errors{2u};
    uint8_t max_errors{2u};
    uint32_t read_length{100u};
    uint32_t number_of_reads{1ULL << 20};
    uint16_t number_of_haplotypes{16u};
};

void run_program(cmd_arguments const & arguments)
{
    std::random_device rd;
    std::mt19937_64 rng(rd());
    std::uniform_int_distribution<uint32_t> read_error_position_dis(0, arguments.read_length - 1);
    std::uniform_int_distribution<uint8_t> dna4_rank_dis(0, 3);

    size_t const number_of_bins{arguments.bin_path.size()};
    uint32_t const reads_per_bin = arguments.number_of_reads / number_of_bins;
    uint32_t const reads_per_haplotype = reads_per_bin / arguments.number_of_haplotypes;
    uint32_t read_counter{};
    // uint64_t bin_counter{};

    std::vector<seqan3::phred42> const quality(arguments.read_length, seqan3::assign_rank_to(40u, seqan3::phred42{}));

    for (auto const & bin_file : arguments.bin_path)
    {
        // std::cerr << "Processing bin " << ++bin_counter << " of " << number_of_bins << '\n';

        // Immediately invoked initialising lambda expession (IIILE).
        std::filesystem::path const out_file = [&]
        {
            std::filesystem::path out_file = arguments.out_path;
            out_file /= bin_file.stem();
            out_file += ".fastq";
            return out_file;
        }();

        seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{bin_file};
        seqan3::sequence_file_output fout{out_file};

        uint16_t haplotype_counter{};
        std::vector<seqan3::dna4> read;

        for (auto const & [seq] : fin)
        {
            uint64_t const reference_length = std::ranges::size(seq);
            std::uniform_int_distribution<uint64_t> read_start_dis(0, reference_length - arguments.read_length);
            for (uint32_t current_read_number = 0; current_read_number < reads_per_haplotype;
                 ++current_read_number, ++read_counter)
            {
                uint64_t const read_start_pos = read_start_dis(rng);
                auto read_slice = seq | seqan3::views::slice(read_start_pos, read_start_pos + arguments.read_length);
                read.assign(read_slice.begin(), read_slice.end());

                for (uint8_t error_count = 0; error_count < arguments.max_errors; ++error_count)
                {
                    uint32_t const error_pos = read_error_position_dis(rng);
                    seqan3::dna4 const current_base = read[error_pos];
                    seqan3::dna4 new_base = current_base;
                    while (new_base == current_base)
                        seqan3::assign_rank_to(dna4_rank_dis(rng), new_base);
                    read[error_pos] = new_base;
                }

                fout.emplace_back(read, std::to_string(read_counter), quality);
            }
            ++haplotype_counter;
        }

        if (haplotype_counter != arguments.number_of_haplotypes)
            std::cerr << "[WARNING] There are not enough / too many haplotypes in the file " << bin_file.string()
                      << '\n'
                      << "[WARNING] Your total read count will be incorrect.\n"
                      << "[WARNING] Haplotypes in file: " << haplotype_counter << '\n'
                      << "[WARNING] Haplotypes expected: " << arguments.number_of_haplotypes << '\n';
    }
}

void initialise_argument_parser(sharg::parser & parser, cmd_arguments & arguments)
{
    parser.info.author = "Enrico Seiler";
    parser.info.author = "enrico.seiler@fu-berlin.de";
    parser.info.short_description = "Generate reads from bins.";
    parser.info.version = "0.0.1";
    parser.info.examples = {"./generate_reads --output ./reads_e2 all_bins.txt"};
    parser.add_positional_option(
        arguments.bin_file_path,
        sharg::config{.description = "Provide a path to a file containing one path to a bin per line.",
                      .validator = sharg::input_file_validator{}});
    parser.add_option(arguments.out_path,
                      sharg::config{.short_id = '\0',
                                    .long_id = "output",
                                    .description = "Provide the base dir where the reads should be written to.",
                                    .required = true,
                                    .validator = sharg::output_directory_validator{}});
    parser.add_option(
        arguments.max_errors,
        sharg::config{.short_id = '\0', .long_id = "max_errors", .description = "The maximum number of errors."});
    parser.add_option(arguments.read_length,
                      sharg::config{.short_id = '\0', .long_id = "read_length", .description = "The read length."});
    parser.add_option(
        arguments.number_of_reads,
        sharg::config{.short_id = '\0', .long_id = "number_of_reads", .description = "The number of reads."});
    parser.add_option(
        arguments.number_of_haplotypes,
        sharg::config{.short_id = '\0', .long_id = "number_of_haplotypes", .description = "The number of haplotypes."});
}

int main(int argc, char ** argv)
{
    sharg::parser myparser{"build_ibf", argc, argv, sharg::update_notifications::off};
    cmd_arguments arguments{};
    initialise_argument_parser(myparser, arguments);
    try
    {
        myparser.parse();
    }
    catch (sharg::parser_error const & ext)
    {
        std::cout << "[Error] " << ext.what() << "\n";
        return -1;
    }

    std::ifstream istrm{arguments.bin_file_path};
    std::string line;
    sharg::input_file_validator validator{};

    while (std::getline(istrm, line))
    {
        if (!line.empty())
        {
            std::filesystem::path bin_path{line};
            validator(bin_path);
            arguments.bin_path.push_back(std::move(bin_path));
        }
    }

    // if (arguments.min_errors > arguments.max_errors)
    //     throw sharg::parser_error{"Minimum number of errors must be at most the maximum number of errors."};

    size_t const number_of_bins{arguments.bin_path.size()};

    if (arguments.max_errors > arguments.read_length)
        throw sharg::parser_error{"Cannot have more errors than the read is long."};

    if (number_of_bins > arguments.number_of_reads)
        throw sharg::parser_error{"Must simulate at least one read per bin."};

    if (arguments.number_of_reads % number_of_bins)
        throw sharg::parser_error{"The number of reads must distribute evenly over the bins."};

    if ((arguments.number_of_reads / number_of_bins) % arguments.number_of_haplotypes)
        throw sharg::parser_error{"The number of reads per bin must evenly distribute over the haplotypes."};

    run_program(arguments);

    return 0;
}

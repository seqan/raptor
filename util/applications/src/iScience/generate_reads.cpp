// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
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

struct config
{
    std::filesystem::path bin_file_path{};
    std::vector<std::filesystem::path> bin_paths{};
    std::filesystem::path out_dir{};
    uint8_t errors{2u};
    uint32_t read_length{100u};
    uint32_t number_of_reads{1ULL << 20};
    uint16_t number_of_haplotypes{16u};
};

void run_program(config const & cfg)
{
    std::random_device rd;
    std::mt19937_64 rng(rd());
    std::uniform_int_distribution<uint32_t> read_error_position_dis(0, cfg.read_length - 1);
    std::uniform_int_distribution<uint8_t> dna4_rank_dis(0, 3);

    size_t const number_of_bins{cfg.bin_paths.size()};
    uint32_t const reads_per_bin = cfg.number_of_reads / number_of_bins;
    uint32_t const reads_per_haplotype = std::max<uint32_t>(1u, reads_per_bin / cfg.number_of_haplotypes);
    uint32_t read_counter{};

    std::vector<seqan3::phred42> const quality(cfg.read_length, seqan3::assign_rank_to(40u, seqan3::phred42{}));
    std::vector<seqan3::dna4> read;

    auto introduce_errors = [&]()
    {
        for (uint8_t error_count{}; error_count < cfg.errors; ++error_count)
        {
            uint32_t const error_pos = read_error_position_dis(rng);
            seqan3::dna4 const current_base = read[error_pos];
            seqan3::dna4 new_base = current_base;
            while (new_base == current_base)
                seqan3::assign_rank_to(dna4_rank_dis(rng), new_base);
            read[error_pos] = new_base;
        }
    };

    for (std::filesystem::path const & bin_file : cfg.bin_paths)
    {
        seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{bin_file};
        seqan3::sequence_file_output fout{(cfg.out_dir / bin_file.filename()).replace_extension(".fastq")};

        uint16_t haplotype_counter{};
        uint32_t reads_per_bin_counter{};

        for (auto && record : fin)
        {
            auto const & seq = record.sequence();
            size_t const reference_length = std::ranges::size(seq);
            std::uniform_int_distribution<size_t> read_start_dis(0, reference_length - cfg.read_length);

            for (uint32_t current_read_number{};
                 current_read_number < reads_per_haplotype && reads_per_bin_counter < reads_per_bin;
                 ++current_read_number, ++read_counter, ++reads_per_bin_counter)
            {
                size_t const read_start_pos = read_start_dis(rng);
                auto read_slice = seqan3::views::slice(seq, read_start_pos, read_start_pos + cfg.read_length);
                read.assign(read_slice.begin(), read_slice.end());

                introduce_errors();

                fout.emplace_back(read, std::to_string(read_counter), quality);
            }
            ++haplotype_counter;
        }
    }
}

void initialise_argument_parser(sharg::parser & parser, config & cfg)
{
    parser.info.author = "Enrico Seiler";
    parser.info.author = "enrico.seiler@fu-berlin.de";
    parser.info.short_description = "Generate reads from bins.";
    parser.info.version = "0.0.1";
    parser.info.examples = {"./generate_reads --output ./reads_e2 all_bins.txt"};
    parser.add_positional_option(cfg.bin_file_path,
                                 sharg::config{.description = "Provide a path to a file containing one path to a "
                                                              "bin per line.",
                                               .validator = sharg::input_file_validator{}});
    parser.add_option(cfg.out_dir,
                      sharg::config{.short_id = '\0',
                                    .long_id = "output",
                                    .description = "Provide the base dir where the reads should be written to.",
                                    .required = true,
                                    .validator = sharg::output_directory_validator{}});
    parser.add_option(
        cfg.errors,
        sharg::config{.short_id = '\0', .long_id = "errors", .description = "The maximum number of errors."});
    parser.add_option(cfg.read_length,
                      sharg::config{.short_id = '\0', .long_id = "read_length", .description = "The read length."});
    parser.add_option(
        cfg.number_of_reads,
        sharg::config{.short_id = '\0', .long_id = "number_of_reads", .description = "The number of reads."});
    parser.add_option(
        cfg.number_of_haplotypes,
        sharg::config{.short_id = '\0', .long_id = "number_of_haplotypes", .description = "The number of haplotypes."});
}

int main(int argc, char ** argv)
{
    sharg::parser myparser{"build_ibf", argc, argv, sharg::update_notifications::off};
    config cfg{};
    initialise_argument_parser(myparser, cfg);
    try
    {
        myparser.parse();
    }
    catch (sharg::parser_error const & ext)
    {
        std::cout << "[Error] " << ext.what() << "\n";
        return -1;
    }

    std::ifstream istrm{cfg.bin_file_path};
    std::string line;
    sharg::input_file_validator validator{};

    while (std::getline(istrm, line))
    {
        if (!line.empty())
        {
            std::filesystem::path bin_path{line};
            validator(bin_path);
            cfg.bin_paths.emplace_back(std::move(bin_path));
        }
    }

    size_t const number_of_bins{cfg.bin_paths.size()};

    if (cfg.errors > cfg.read_length)
        throw sharg::parser_error{"Cannot have more errors than the read is long."};

    if (number_of_bins > cfg.number_of_reads)
        throw sharg::parser_error{"Must simulate at least one read per bin."};

    if (cfg.number_of_reads % number_of_bins)
        throw sharg::parser_error{"The number of reads must distribute evenly over the bins."};

    run_program(cfg);

    return 0;
}

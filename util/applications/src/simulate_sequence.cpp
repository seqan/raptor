// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <random>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sequence_file/output.hpp>

#include <raptor/argument_parsing/validators.hpp>

struct config
{
    size_t sequence_length{};
    size_t seed{};
    std::filesystem::path output{"sequence.fa"};
};

inline void set_up_parser(sharg::parser & parser, config & cfg)
{
    parser.info.version = "1.0.0";
    parser.info.author = "Enrico Seiler";
    parser.info.email = "enrico.seiler@fu-berlin.de";
    parser.info.short_description = "Simulates a random DNA sequence of a given length.";

    parser.add_subsection("Main options:");
    parser.add_option(cfg.sequence_length,
                      sharg::config{.short_id = '\0',
                                    .long_id = "length",
                                    .description = "The length of the sequence to simulate. ",
                                    .required = true});
    parser.add_option(cfg.seed,
                      sharg::config{
                          .short_id = '\0',
                          .long_id = "seed",
                          .description = "The seed for the random number generator. ",
                      });
    parser.add_option(cfg.output,
                      sharg::config{.short_id = '\0',
                                    .long_id = "output",
                                    .description = "The output. ",
                                    .validator = sharg::output_file_validator{{"fa", "fasta"}}});
}

inline auto simulate_sequence(config const & cfg)
{
    std::mt19937_64 rng{cfg.seed};
    std::uniform_int_distribution<size_t> dna4_rank_dis{0u, 3u};

    std::vector<seqan3::dna4> sequence(cfg.sequence_length);
    for (auto & letter : sequence)
        letter.assign_rank(dna4_rank_dis(rng));

    return sequence;
}

int main(int argc, char const * argv[])
{
    sharg::parser parser{"simulate_sequence", argc, argv, sharg::update_notifications::off};
    parser.info.version = "1.0.0";

    config cfg{};
    set_up_parser(parser, cfg);

    try
    {
        parser.parse();
    }
    catch (sharg::parser_error const & ext)
    {
        std::cerr << "[ERROR] " << ext.what() << '\n';
        return -1;
    }

    seqan3::sequence_file_output fout{cfg.output};
    fout.emplace_back(simulate_sequence(cfg), std::string_view{"sequence"});
}

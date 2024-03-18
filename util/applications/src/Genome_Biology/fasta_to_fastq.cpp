// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <sharg/all.hpp>

#include <seqan3/io/sequence_file/all.hpp>

struct dna4_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

int main(int argc, char ** argv)
{
    std::filesystem::path input_file{};
    std::filesystem::path output_file{};

    sharg::parser parser{"fasta_to_fastq", argc, argv, sharg::update_notifications::off};
    parser.add_option(
        input_file,
        sharg::config{.short_id = '\0', .long_id = "input", .description = "Input FASTA file.", .required = true});
    parser.add_option(
        output_file,
        sharg::config{.short_id = '\0', .long_id = "output", .description = "OUTPUT FASTQ file.", .required = true});
    parser.parse();

    seqan3::sequence_file_input<dna4_traits> fin{input_file};
    seqan3::sequence_file_output fout{output_file};
    std::vector<seqan3::phred42> quality{};
    seqan3::phred42 const default_quality = seqan3::assign_rank_to(40u, seqan3::phred42{});

    for (auto && record : fin)
    {
        quality.resize(std::ranges::size(record.sequence()), default_quality);
        fout.emplace_back(record.sequence(), record.id(), quality);
    }
}

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>

struct dna4_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

int main(int argc, char ** argv)
{
    std::filesystem::path input_file{};
    std::filesystem::path output_file{};

    seqan3::argument_parser parser{"fasta_to_fastq", argc, argv, seqan3::update_notifications::off};
    parser.add_option(input_file, '\0', "input", "Input FASTA file.", seqan3::option_spec::required);
    parser.add_option(output_file, '\0', "output", "OUTPUT FASTQ file.", seqan3::option_spec::required);
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

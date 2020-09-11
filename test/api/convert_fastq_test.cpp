#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "fastq_conversion.hpp"

TEST(group1, out_empty)
{
    std::string expected{"> seq1\nACGTTTGATTCGCG\n> seq2\nTCGGGGGATTCGCG\n"};
    testing::internal::CaptureStdout();
    convert_fastq(DATADIR"in.fastq", "");
    std::string std_cout = testing::internal::GetCapturedStdout();
    EXPECT_EQ(expected, std_cout);
}

TEST(group1, out_not_empty)
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory
    convert_fastq(DATADIR"in.fastq", tmp_dir/"out.fasta");                 // create out.fasta

    // Check if out.fasta is correct
    using seqan3::operator""_dna5;
    std::vector<seqan3::dna5_vector> expected_seqs = {"ACGTTTGATTCGCG"_dna5, "TCGGGGGATTCGCG"_dna5};
    std::vector<std::string> expected_ids{"seq1", "seq2"};
    std::vector<seqan3::dna5_vector> seqs{};
    std::vector<std::string> ids{};
    seqan3::sequence_file_input fin{tmp_dir/"out.fasta"};

    for (auto & [seq, id, qual] : fin)
    {
        ids.push_back(id);
        seqs.push_back(seq);
    }
    EXPECT_EQ(seqs, expected_seqs);
    EXPECT_EQ(ids, expected_ids);
}

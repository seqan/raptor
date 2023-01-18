#include <seqan3/core/debug_stream.hpp>

#include <raptor/hierarchical_interleaved_bloom_filter.hpp>

// For this example we have two input fasta files with the following three sequences (= user bins):
// example1.fasta:
// ```fasta
// >chr1
// AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGCGTTCATTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
// ```
// example2.fasta:
// ```fasta
// >chr2
// AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGCGTCATTAA
// >chr3
// AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
// ```

// 1-level-IBF: |chr1|chr2,chr3|
// 2-level-IBF:      |chr2|chr3|

int main()
{
    raptor::hierarchical_interleaved_bloom_filter hibf{};

    hibf.user_bins.set_user_bin_count(3);

    hibf.user_bins.filename_of_user_bin(0) = "path/example1.fasta";
    hibf.user_bins.filename_of_user_bin(1) = "path/example2.fasta";
    hibf.user_bins.filename_of_user_bin(2) = "path/example2.fasta";

    seqan3::debug_stream << "Filename of user bin 0: " << hibf.user_bins.filename_of_user_bin(0) << '\n';
    seqan3::debug_stream << "Filename of user bin 1: " << hibf.user_bins.filename_of_user_bin(1) << '\n';
    seqan3::debug_stream << "Filename of user bin 2: " << hibf.user_bins.filename_of_user_bin(2) << '\n';
}

// Prints out:
// Filename of user bin 0: path/example1.fasta
// Filename of user bin 1: path/example2.fasta
// Filename of user bin 2: path/example2.fasta

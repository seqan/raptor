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

    hibf.user_bins.set_ibf_count(2);

    hibf.user_bins.bin_indices_of_ibf(0) = {1, 2, 3};
    hibf.user_bins.bin_indices_of_ibf(1) = {2, 3};

    seqan3::debug_stream << "User bin indices of 1-level-IBF: " << hibf.user_bins.bin_indices_of_ibf(0) << '\n';
    seqan3::debug_stream << "User bin indices of 2-level-IBF: " << hibf.user_bins.bin_indices_of_ibf(1) << '\n';
}

// Prints out:
// User bin indices of 1-level-IBF: [1,2,3]
// User bin indices of 2-level-IBF: [2,3]

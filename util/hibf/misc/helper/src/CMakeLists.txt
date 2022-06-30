cmake_minimum_required (VERSION 3.8)

add_executable ("fasta_to_fastq" fasta_to_fastq.cpp)
target_link_libraries ("fasta_to_fastq" PUBLIC seqan3::seqan3)

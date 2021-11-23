cmake_minimum_required (VERSION 3.8)

add_executable ("${PROJECT_NAME}" fasta_to_fastq.cpp)
target_link_libraries ("${PROJECT_NAME}" PUBLIC seqan3::seqan3)


cmake_minimum_required (VERSION 3.18)

if (NOT TARGET fasta_to_fastq)
    add_executable ("fasta_to_fastq" fasta_to_fastq.cpp)
    target_link_libraries ("fasta_to_fastq" "utility_common")
endif ()

if (NOT TARGET normalise_yara_truth_file)
    add_executable ("normalise_yara_truth_file" normalise_yara_truth_file.cpp)
    target_link_libraries ("normalise_yara_truth_file" "utility_common")
endif ()

if (NOT TARGET compare_to_truth)
    add_executable ("compare_to_truth" compare_to_truth.cpp)
    target_link_libraries ("compare_to_truth" "utility_common")
endif ()

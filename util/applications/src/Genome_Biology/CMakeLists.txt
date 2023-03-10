cmake_minimum_required (VERSION 3.18)

add_executable ("normalise_mantis_output" normalise_mantis_output.cpp)
target_link_libraries ("normalise_mantis_output" "utility_common")

add_executable ("compare_mantis_raptor_output" compare_mantis_raptor_output.cpp)
target_link_libraries ("compare_mantis_raptor_output" "utility_common")

add_executable ("normalise_raptor_output" normalise_raptor_output.cpp)
target_link_libraries ("normalise_raptor_output" "utility_common")

add_executable ("fasta_to_fastq" fasta_to_fastq.cpp)
target_link_libraries ("fasta_to_fastq" "utility_common")

add_executable ("normalise_yara_truth_file" normalise_yara_truth_file.cpp)
target_link_libraries ("normalise_yara_truth_file" "utility_common")

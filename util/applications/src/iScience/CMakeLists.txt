cmake_minimum_required (VERSION 3.18)

add_executable ("generate_reads" generate_reads.cpp)
target_link_libraries ("generate_reads" "utility_common")

add_executable ("split_sequence" split_sequence.cpp)
target_link_libraries ("split_sequence" "utility_common")

add_executable ("count_minimiser" ${RAPTOR_CLONE_DIR}/src/argument_parsing/parse_bin_path.cpp count_minimiser.cpp)
target_link_libraries ("count_minimiser" "utility_common")

add_executable ("threshold_info"
                ${RAPTOR_CLONE_DIR}/src/threshold/multiple_error_model.cpp
                ${RAPTOR_CLONE_DIR}/src/threshold/one_error_model.cpp
                ${RAPTOR_CLONE_DIR}/src/threshold/one_indirect_error_model.cpp
                ${RAPTOR_CLONE_DIR}/src/threshold/pascal_row.cpp
                ${RAPTOR_CLONE_DIR}/src/threshold/precompute_correction.cpp
                ${RAPTOR_CLONE_DIR}/src/threshold/precompute_threshold.cpp
                threshold_info.cpp
)
target_link_libraries ("threshold_info" "utility_common")

add_executable ("apply_taxsbp" apply_taxsbp.cpp)
target_link_libraries ("apply_taxsbp" "utility_common")

add_executable ("generate_reads_refseq" generate_reads_refseq.cpp)
target_link_libraries ("generate_reads_refseq" "utility_common")

add_executable ("ibf_fpr" ibf_fpr.cpp)
target_link_libraries ("ibf_fpr" "utility_common")

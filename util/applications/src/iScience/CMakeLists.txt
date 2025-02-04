# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25...3.30)

if (NOT TARGET split_sequence)
    add_executable ("split_sequence" split_sequence.cpp)
    target_link_libraries ("split_sequence" "utility_common")
endif ()

if (NOT TARGET count_minimiser)
    add_executable ("count_minimiser" ${Raptor_SOURCE_DIR}/src/argument_parsing/parse_bin_path.cpp count_minimiser.cpp)
    target_link_libraries ("count_minimiser" "utility_common" chopper::layout)
endif ()

if (NOT TARGET threshold_info)
    add_executable ("threshold_info"
                    ${Raptor_SOURCE_DIR}/src/threshold/multiple_error_model.cpp
                    ${Raptor_SOURCE_DIR}/src/threshold/one_error_model.cpp
                    ${Raptor_SOURCE_DIR}/src/threshold/one_indirect_error_model.cpp
                    ${Raptor_SOURCE_DIR}/src/threshold/pascal_row.cpp
                    ${Raptor_SOURCE_DIR}/src/threshold/precompute_correction.cpp
                    ${Raptor_SOURCE_DIR}/src/threshold/precompute_threshold.cpp
                    threshold_info.cpp
    )
    target_link_libraries ("threshold_info" "utility_common")
endif ()

if (NOT TARGET apply_taxsbp)
    add_executable ("apply_taxsbp" apply_taxsbp.cpp)
    target_link_libraries ("apply_taxsbp" "utility_common")
endif ()

if (NOT TARGET generate_reads)
    add_executable ("generate_reads"
                    generate_reads.cpp
                    generate_reads/apply_weights.cpp
                    generate_reads/infer_weights_from_file_size.cpp
                    generate_reads/infer_weights_from_hll_sketches.cpp
                    generate_reads/infer_weights_from_kmer_counts.cpp
                    generate_reads/infer_weights_from_uniform_distribution.cpp
                    generate_reads/initialise_argument_parser.cpp
                    generate_reads/parse_input.cpp
                    generate_reads/print_weights.cpp
                    generate_reads/simulate_reads.cpp
    )
    target_link_libraries ("generate_reads" "utility_common")
endif ()

if (NOT TARGET ibf_fpr)
    add_executable ("ibf_fpr" ibf_fpr.cpp)
    target_link_libraries ("ibf_fpr" "utility_common")
endif ()

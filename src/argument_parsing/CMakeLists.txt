# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25...3.30)

if (TARGET raptor::argument_parsing)
    return ()
endif ()

add_library ("raptor_argument_parsing" STATIC
             build_arguments.cpp
             build_parsing.cpp
             compute_bin_size.cpp
             parse_bin_path.cpp
             prepare_parsing.cpp
             search_arguments.cpp
             search_parsing.cpp
             update_parsing.cpp
             upgrade_parsing.cpp
)

set_source_files_properties (search_parsing.cpp
                             PROPERTIES COMPILE_DEFINITIONS
                                        "FPGA_WINDOWS=${WINDOW_SIZE_STRING};FPGA_KMERS=${MIN_IBF_K_STRING};FPGA_BINS=${BIN_COUNT_STRING};FPGA_KERNELS=${KERNEL_COPYS_STRING}"
)

target_link_libraries ("raptor_argument_parsing" PUBLIC "raptor::interface")
add_library (raptor::argument_parsing ALIAS raptor_argument_parsing)

# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25...3.30)

if (TARGET raptor::search)
    return ()
endif ()

add_library ("raptor_search" STATIC raptor_search.cpp search_hibf.cpp search_ibf.cpp search_partitioned_ibf.cpp)
target_link_libraries ("raptor_search" PUBLIC "raptor::interface")

if (RAPTOR_FPGA)
    add_subdirectory (fpga)
    target_link_libraries ("raptor_search" PUBLIC "raptor::search::fpga")
endif ()

add_library (raptor::search ALIAS raptor_search)

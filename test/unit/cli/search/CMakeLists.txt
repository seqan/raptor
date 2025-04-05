# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25...3.30)

raptor_add_unit_test (search_hibf_test.cpp)
raptor_add_unit_test (search_hibf_preprocessing_test.cpp)
raptor_add_unit_test (search_ibf_preprocessing_test.cpp)
raptor_add_unit_test (search_ibf_test.cpp)

if (RAPTOR_FPGA)
    raptor_add_unit_test (search_ibf_fpga_test.cpp)
endif ()

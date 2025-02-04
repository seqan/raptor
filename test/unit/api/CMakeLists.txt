# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25...3.30)

raptor_add_unit_test (compute_bin_size.cpp)
raptor_add_unit_test (issue_142.cpp)
raptor_add_unit_test (formatted_bytes.cpp)
raptor_add_unit_test (index_size.cpp)
raptor_add_unit_test (memory_usage.cpp)
raptor_add_unit_test (threshold.cpp)
raptor_add_unit_test (validate_shape.cpp)

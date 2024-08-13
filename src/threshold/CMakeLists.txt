# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25)

if (TARGET raptor::threshold)
    return ()
endif ()

add_library ("raptor_threshold" STATIC
             multiple_error_model.cpp
             one_error_model.cpp
             one_indirect_error_model.cpp
             pascal_row.cpp
             precompute_correction.cpp
             precompute_threshold.cpp
             threshold.cpp
)

target_link_libraries ("raptor_threshold" PUBLIC "raptor::interface")
add_library (raptor::threshold ALIAS raptor_threshold)

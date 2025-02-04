# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.15...3.30)

if (TARGET raptor::layout)
    return ()
endif ()

add_library ("raptor_layout" STATIC raptor_layout.cpp)
target_link_libraries ("raptor_layout" PUBLIC "raptor::interface" "chopper::chopper")
add_library (raptor::layout ALIAS raptor_layout)

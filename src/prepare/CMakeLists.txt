# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25)

if (TARGET raptor_prepare)
    return ()
endif ()

add_library ("raptor_prepare" STATIC compute_minimiser.cpp)

target_link_libraries ("raptor_prepare" PUBLIC "raptor::raptor")

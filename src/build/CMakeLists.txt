# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25)

if (TARGET raptor_build)
    return ()
endif ()

add_library ("raptor_build" STATIC build_hibf.cpp build_ibf.cpp max_count_per_partition.cpp raptor_build.cpp)

target_link_libraries ("raptor_build" PUBLIC "raptor::raptor" "raptor_prepare" "seqan::hibf")

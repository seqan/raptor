# SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25...3.30)

if (TARGET raptor::inspect)
    return ()
endif ()

add_library ("raptor_inspect" STATIC index_metadata.cpp json_writer.cpp raptor_inspect.cpp write_metadata.cpp)
target_link_libraries ("raptor_inspect" PUBLIC "raptor::interface")
add_library (raptor::inspect ALIAS raptor_inspect)

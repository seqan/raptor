# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25)

if (TARGET raptor::update)
    return ()
endif ()

add_library ("raptor_update" STATIC delete_user_bins.cpp dump_index.cpp insert/get_location.cpp
                                    insert/insert_tb_and_parents.cpp insert_user_bin.cpp raptor_update.cpp
)
target_link_libraries ("raptor_update" PUBLIC "raptor::interface")
add_library (raptor::update ALIAS raptor_update)

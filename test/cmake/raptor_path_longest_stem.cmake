# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25...3.30)

# A compatible function for cmake < 3.20 that basically returns `cmake_path (GET <filename> STEM LAST_ONLY <out_var>)`
function (raptor_path_longest_stem out_var filename)
    cmake_path (GET filename STEM LAST_ONLY result)

    set ("${out_var}"
         "${result}"
         PARENT_SCOPE
    )
endfunction ()

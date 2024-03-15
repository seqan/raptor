# --------------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
# --------------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.25)

# A compatible function for cmake < 3.20 that basically returns `cmake_path (GET <filename> STEM LAST_ONLY <out_var>)`
function (raptor_path_longest_stem out_var filename)
    cmake_path (GET filename STEM LAST_ONLY result)

    set ("${out_var}"
         "${result}"
         PARENT_SCOPE
    )
endfunction ()

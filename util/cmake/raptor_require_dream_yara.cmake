# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

macro (raptor_require_dream_yara)
    set (dream_yara_args ${SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS})

    include (ExternalProject)
    ExternalProject_Add (
        dream_yara
        PREFIX dream_yara
        GIT_REPOSITORY "https://github.com/seqan/dream_yara.git"
        GIT_TAG "raptor_ibf"
        GIT_SHALLOW true
        GIT_SUBMODULES_RECURSE true
        SOURCE_DIR "${PROJECT_BINARY_DIR}/dream_yara/dream_yara"
        CMAKE_ARGS "${dream_yara_args}"
        INSTALL_COMMAND ""
    )
    unset (dream_yara_args)

    ExternalProject_Get_property (dream_yara BINARY_DIR)
    install (FILES "${BINARY_DIR}/bin/dream_yara_build_filter"
             DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
             PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                         GROUP_READ             GROUP_EXECUTE
                         WORLD_READ             WORLD_EXECUTE)
    install (FILES "${BINARY_DIR}/bin/dream_yara_indexer"
             DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
             PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                         GROUP_READ             GROUP_EXECUTE
                         WORLD_READ             WORLD_EXECUTE)
    install (FILES "${BINARY_DIR}/bin/dream_yara_mapper"
             DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
             PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                         GROUP_READ             GROUP_EXECUTE
                         WORLD_READ             WORLD_EXECUTE)
endmacro ()

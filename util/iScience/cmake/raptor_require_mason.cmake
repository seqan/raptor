# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

macro (raptor_require_mason)
    block (SCOPE_FOR VARIABLES)
    set (mason_args ${SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS})
    list (APPEND mason_args "-DSEQAN_BUILD_SYSTEM=APP:mason2")

    include (ExternalProject)
    ExternalProject_Add (mason
                         PREFIX ${CMAKE_BINARY_DIR}/mason
                         GIT_REPOSITORY "https://github.com/seqan/seqan.git"
                         GIT_TAG "main"
                         GIT_SHALLOW true
                         SOURCE_DIR "${PROJECT_BINARY_DIR}/mason/seqan"
                         CMAKE_ARGS "${mason_args}"
                         BUILD_COMMAND ${CMAKE_COMMAND} --build . --target mason_variator
                         INSTALL_COMMAND ""
    )
    unset (mason_args)

    ExternalProject_Get_Property (mason BINARY_DIR)

    install (FILES "${BINARY_DIR}/bin/mason_variator"
             DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
             PERMISSIONS OWNER_READ
                         OWNER_WRITE
                         OWNER_EXECUTE
                         GROUP_READ
                         GROUP_EXECUTE
                         WORLD_READ
                         WORLD_EXECUTE
    )

    add_custom_command (TARGET mason
                        POST_BUILD
                        COMMAND ${CMAKE_COMMAND} -E copy ${BINARY_DIR}/bin/mason_variator
                                ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/mason_variator
    )
    endblock ()
endmacro ()

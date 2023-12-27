# --------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
# --------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.21)

find_path (RAPTOR_MODULE_PATH "configuration.cmake" HINTS "${CMAKE_CURRENT_LIST_DIR}/cmake/")
list (APPEND CMAKE_MODULE_PATH "${RAPTOR_MODULE_PATH}")

include (version)

project (Raptor
         LANGUAGES CXX
         VERSION "${RAPTOR_PROJECT_VERSION}"
         DESCRIPTION "Raptor -- A fast and space-efficient pre-filter for querying very large collections of nucleotide sequences"
         HOMEPAGE_URL "https://github.com/seqan/raptor"
)

if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE
         Release
         CACHE STRING "Choose the type of build, options are: Debug Coverage Release RelWithDebInfo MinSizeRel." FORCE
    )
endif ()

include (configuration)

add_subdirectory (src)

option (INSTALL_RAPTOR "Enable installation of Raptor." ON)

if (INSTALL_RAPTOR)
    include (GNUInstallDirs)
    install (TARGETS raptor RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
endif ()

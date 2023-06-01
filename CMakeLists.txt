# --------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
# --------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.18)

find_path (RAPTOR_MODULE_PATH "raptor-config.cmake" HINTS "${CMAKE_CURRENT_LIST_DIR}/build_system/")
list (APPEND CMAKE_MODULE_PATH "${RAPTOR_MODULE_PATH}")

include (raptor-config-version)

project (raptor
         LANGUAGES CXX
         VERSION "${RAPTOR_PROJECT_VERSION}"
         DESCRIPTION "Raptor -- A fast and space-efficient pre-filter for querying very large collections of nucleotide sequences"
         HOMEPAGE_URL "https://github.com/seqan/raptor"
)

#test
find_package (Raptor REQUIRED HINTS ${RAPTOR_MODULE_PATH})

option (BUILD_RAPTOR "Enable building of Raptor." ON)

if (BUILD_RAPTOR)
    add_subdirectory (src)
endif ()

option (INSTALL_RAPTOR_HEADERS "Enable installation of Raptor headers." OFF)

if (INSTALL_RAPTOR_HEADERS)
    include (raptor-install)
    include (raptor-package)
endif ()

# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

file (STRINGS "${CMAKE_CURRENT_LIST_DIR}/../include/raptor/version.hpp" RAPTOR_VERSION_HPP
      REGEX "#define RAPTOR_VERSION_(MAJOR|MINOR|PATCH)"
)
string (REGEX REPLACE "#define RAPTOR_VERSION_(MAJOR|MINOR|PATCH) " "" RAPTOR_VERSION "${RAPTOR_VERSION_HPP}")
string (REGEX REPLACE ";" "." RAPTOR_VERSION "${RAPTOR_VERSION}")

file (STRINGS "${CMAKE_CURRENT_LIST_DIR}/../include/raptor/version.hpp" RAPTOR_RELEASE_CANDIDATE_HPP
      REGEX "#define RAPTOR_RELEASE_CANDIDATE "
)
string (REGEX REPLACE "#define RAPTOR_RELEASE_CANDIDATE " "" RAPTOR_RELEASE_CANDIDATE_VERSION
                      "${RAPTOR_RELEASE_CANDIDATE_HPP}"
)

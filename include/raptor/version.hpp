// --------------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides RAPTOR version macros and global variables.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <cstdint>

//!\brief The major version as MACRO.
#define RAPTOR_VERSION_MAJOR 3
//!\brief The minor version as MACRO.
#define RAPTOR_VERSION_MINOR 1
//!\brief The patch version as MACRO.
#define RAPTOR_VERSION_PATCH 0
//!\brief The release candidate number. 0 means stable release, >= 1 means release candidate.
#define RAPTOR_RELEASE_CANDIDATE 1

//!\brief The full version as MACRO (number).
#define RAPTOR_VERSION (RAPTOR_VERSION_MAJOR * 10000 + RAPTOR_VERSION_MINOR * 100 + RAPTOR_VERSION_PATCH)

/*!\brief Converts a number to a string. Preprocessor needs this indirection to
 * properly expand the values to strings.
 */
#define RAPTOR_VERSION_CSTRING_HELPER_STR(str) #str

//!\brief Converts version numbers to string.
#define RAPTOR_VERSION_CSTRING_HELPER_FUNC(MAJOR, MINOR, PATCH)                                                        \
    RAPTOR_VERSION_CSTRING_HELPER_STR(MAJOR)                                                                           \
    "." RAPTOR_VERSION_CSTRING_HELPER_STR(MINOR) "." RAPTOR_VERSION_CSTRING_HELPER_STR(PATCH)

#if (RAPTOR_RELEASE_CANDIDATE > 0)
//!\brief A helper function that expands to a suitable release candidate suffix.
#    define RAPTOR_RELEASE_CANDIDATE_HELPER(RC) "-rc." RAPTOR_VERSION_CSTRING_HELPER_STR(RC)
#else
//!\brief A helper function that expands to a suitable release candidate suffix.
#    define RAPTOR_RELEASE_CANDIDATE_HELPER(RC) ""
#endif

//!\brief The full version as null terminated string.
#define RAPTOR_VERSION_CSTRING                                                                                         \
    RAPTOR_VERSION_CSTRING_HELPER_FUNC(RAPTOR_VERSION_MAJOR, RAPTOR_VERSION_MINOR, RAPTOR_VERSION_PATCH)               \
    RAPTOR_RELEASE_CANDIDATE_HELPER(RAPTOR_RELEASE_CANDIDATE)

namespace raptor
{

//!\brief The major version.
constexpr uint8_t raptor_version_major = RAPTOR_VERSION_MAJOR;
//!\brief The minor version.
constexpr uint8_t raptor_version_minor = RAPTOR_VERSION_MINOR;
//!\brief The patch version.
constexpr uint8_t raptor_version_patch = RAPTOR_VERSION_PATCH;

//!\brief The full version as `std::size_t`.
constexpr std::size_t raptor_version = RAPTOR_VERSION;

//!\brief The full version as null terminated string.
constexpr char const * raptor_version_cstring = RAPTOR_VERSION_CSTRING;

} // namespace raptor

#undef RAPTOR_VERSION_CSTRING_HELPER_STR
#undef RAPTOR_VERSION_CSTRING_HELPER_FUNC
#undef RAPTOR_RELEASE_CANDIDATE_HELPER

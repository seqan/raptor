// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::formatted_peak_ram.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <cstddef> // for size_t
#include <string>  // for basic_string, string

#include <raptor/argument_parsing/formatted_bytes.hpp> // for formatted_bytes

#if __has_include(<sys/resource.h>)
#    include <sys/resource.h> // for rusage, getrusage, RUSAGE_SELF
#endif

namespace raptor
{

#if __has_include(<sys/resource.h>)
// Returns -1 if not available. Actually returns bytes instead of KiB on macOS.
inline long peak_ram_in_KiB()
{
    rusage usage;
#    if __APPLE__
    return getrusage(RUSAGE_SELF, &usage) == 0 ? (usage.ru_maxrss >> 10) : -1L;
#    else
    return getrusage(RUSAGE_SELF, &usage) == 0 ? usage.ru_maxrss : -1L;
#    endif
}
#else
inline long peak_ram_in_KiB()
{
    return -1L;
}
#endif

[[nodiscard]] inline std::string formatted_peak_ram()
{
    long const peak_ram_KiB = peak_ram_in_KiB();
    if (peak_ram_KiB == -1L)
        return {": Not available"}; // GCOVR_EXCL_LINE

    return formatted_bytes(static_cast<size_t>(peak_ram_KiB) << 10);
}

} // namespace raptor

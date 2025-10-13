// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::get_cpu_time.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <algorithm> // for min
#include <cassert>   // for assert
#include <chrono>    // for duration, operator+, microseconds, seconds
#include <cstdint>   // for uint8_t

#if __has_include(<sys/resource.h>)
#    include <sys/resource.h> // for rusage, timeval, getrusage, RUSAGE_SELF
#endif

namespace raptor
{

struct cpu_time_t
{
    double user_time_in_seconds{};
    double system_time_in_seconds{};

    void operator-=(double const diff) noexcept
    {
        user_time_in_seconds -= diff;
    }

    [[nodiscard]] double cpu_usage_in_percent(double const elapsed_time_in_seconds,
                                              uint8_t const threads) const noexcept
    {
        assert(is_valid());
        assert(elapsed_time_in_seconds > 0.0);
        return std::min<double>(100.0 * threads,
                                100.0 * (user_time_in_seconds + system_time_in_seconds) / elapsed_time_in_seconds);
    }

    [[nodiscard]] bool is_valid() const noexcept
    {
        return (user_time_in_seconds + system_time_in_seconds) >= 0.0;
    }
};

#if __has_include(<sys/resource.h>)
// Returns cpu_time_t{-1.0, -1.0} if not available.
[[nodiscard]] inline cpu_time_t get_cpu_time()
{
    rusage usage;

    if (getrusage(RUSAGE_SELF, &usage) != 0)
        return cpu_time_t{-1.0, -1.0}; // GCOVR_EXCL_LINE

    std::chrono::duration user_time =
        std::chrono::seconds{usage.ru_utime.tv_sec} + std::chrono::microseconds{usage.ru_utime.tv_usec};
    std::chrono::duration system_time =
        std::chrono::seconds{usage.ru_stime.tv_sec} + std::chrono::microseconds{usage.ru_stime.tv_usec};

    return cpu_time_t{.user_time_in_seconds = std::chrono::duration<double>(user_time).count(),
                      .system_time_in_seconds = std::chrono::duration<double>(system_time).count()};
}
#else
[[nodiscard]] inline cpu_time_t get_cpu_time()
{
    return cpu_time_t{-1.0, -1.0};
}
#endif

} // namespace raptor

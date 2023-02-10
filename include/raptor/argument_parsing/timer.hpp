// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides raptor::timer.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <atomic>
#include <cassert>
#include <chrono>

namespace raptor
{

enum class concurrent
{
    no,
    yes
};

template <concurrent concurrency>
class timer
{
private:
    static constexpr bool is_concurrent{concurrency == concurrent::yes};

    using rep_t =
        std::conditional_t<is_concurrent, std::atomic<std::chrono::steady_clock::rep>, std::chrono::steady_clock::rep>;

    template <concurrent concurrency_>
    friend class timer;

public:
    timer() = default;
    timer(timer &&) = default;
    timer & operator=(timer &&) = default;
    ~timer() = default;

    timer(timer const & other)
        requires (!is_concurrent)
    = default;
    timer & operator=(timer const & other)
        requires (!is_concurrent)
    = default;

    timer(timer const & other)
        requires is_concurrent
        : start_{other.start_}, stop_{other.stop_}, ticks{other.ticks.load()}
    {}
    timer & operator=(timer const & other)
        requires is_concurrent
    {
        start_ = other.start_;
        stop_ = other.stop_;
        ticks = other.ticks.load();
        return *this;
    }

    void start()
    {
        start_ = std::chrono::steady_clock::now();
    }

    void stop()
    {
        stop_ = std::chrono::steady_clock::now();
        assert(stop_ >= start_);
        ticks += (stop_ - start_).count();
    }

    template <concurrent concurrency_>
    void operator+=(timer<concurrency_> const & other)
    {
        ticks += other.ticks;
    }

    double in_seconds() const
        requires is_concurrent
    {
        return std::chrono::duration<double>(std::chrono::steady_clock::duration{ticks.load()}).count();
    }

    double in_seconds() const
        requires (!is_concurrent)
    {
        return std::chrono::duration<double>(std::chrono::steady_clock::duration{ticks}).count();
    }

private:
    std::chrono::steady_clock::time_point start_{std::chrono::time_point<std::chrono::steady_clock>::max()};
    std::chrono::steady_clock::time_point stop_{};
    rep_t ticks{};
};

} // namespace raptor

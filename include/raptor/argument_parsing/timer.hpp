// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <atomic>
#include <cassert>
#include <chrono>

namespace raptor
{

class timer
{
public:
    timer() = default;
    timer(timer const & other) : start_{other.start_}, stop_{other.stop_} // this will have ticks = 0
    {}
    timer(timer &&) = default;
    timer & operator=(timer const & other) // this will have ticks = 0
    {
        start_ = other.start_;
        stop_ = other.stop_;
        return *this;
    }
    timer & operator=(timer &&) = default;
    ~timer() = default;

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

    void operator+=(timer const & other)
    {
        ticks += other.ticks;
    }

    double in_seconds() const
    {
        return std::chrono::duration<double>(std::chrono::steady_clock::duration{ticks.load()}).count();
    }

private:
    std::chrono::steady_clock::time_point start_{std::chrono::time_point<std::chrono::steady_clock>::max()};
    std::chrono::steady_clock::time_point stop_{};
    std::atomic<std::chrono::steady_clock::rep> ticks{};
};

} // namespace raptor

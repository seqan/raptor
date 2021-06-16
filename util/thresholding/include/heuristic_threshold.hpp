// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/eseiler/minimizer_thresholds/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides heuristic_threshold.
 */

#pragma once

#include <algorithm>
#include <iostream>
#include <stdint.h>
#include <vector>

#include "minimizer.hpp"

struct heuristic_threshold
{
    std::vector<uint64_t> minimizer_begin;
    std::vector<uint64_t> minimizer_end;
    std::vector<uint32_t> coverage;
    std::vector<uint64_t> coverage_begin;
    std::vector<uint64_t> coverage_end;

    heuristic_threshold() = default;
    heuristic_threshold(heuristic_threshold const &) = default;
    heuristic_threshold & operator=(heuristic_threshold const &) = default;
    heuristic_threshold(heuristic_threshold &&) = default;
    heuristic_threshold & operator=(heuristic_threshold &&) = default;
    ~heuristic_threshold() = default;

    heuristic_threshold(minimizer const & mini) :
        minimizer_begin{mini.minimizer_begin}, minimizer_end{mini.minimizer_end}
    {}

    inline void compute_coverage()
    {
        uint64_t begin_pos{1};
        uint64_t end_pos{0};

        auto unique_minimizer_end = minimizer_end;
        unique_minimizer_end.erase(unique(unique_minimizer_end.begin(), unique_minimizer_end.end()), unique_minimizer_end.end());

        auto unique_minimizer_begin = minimizer_begin;
        unique_minimizer_begin.erase(unique(unique_minimizer_begin.begin(), unique_minimizer_begin.end()), unique_minimizer_begin.end());

        coverage_begin.push_back(unique_minimizer_begin[0]);
        coverage.push_back(1);

        while ((begin_pos < unique_minimizer_begin.size()) || (end_pos < unique_minimizer_end.size()))
        {
            uint64_t begin = begin_pos < unique_minimizer_begin.size() ? unique_minimizer_begin[begin_pos] : 0xFFFFFFFFFFFFFFFFULL;
            uint64_t end = unique_minimizer_end[end_pos];
            // Overlap
            if (begin < end)
            {
                coverage_end.push_back(begin - 1);
                coverage_begin.push_back(begin);
                coverage.push_back(coverage.back() + 1);
                ++begin_pos;
            }
            // Flatten consecutive positions, where one kmer ends and other one starts
            if (begin == end)
            {
                coverage_end.push_back(begin - 1);
                coverage_begin.push_back(begin);
                coverage.push_back(coverage.back() + 1);
                while (unique_minimizer_begin[begin_pos] == unique_minimizer_end[end_pos])
                {
                    ++begin_pos;
                    ++end_pos;
                }
                --end_pos;
            }
            // Kmer ends
            if (end < begin)
            {
                coverage_end.push_back(end);
                coverage_begin.push_back(end + 1);
                coverage.push_back(coverage.back() - 1);
                ++end_pos;
            }
        }
        coverage_begin.pop_back();
        coverage.pop_back();
    }

    // t = text length, e = errors
    inline uint32_t threshold(uint16_t errors)
    {
        uint32_t destroyed{0};
        uint32_t available{static_cast<uint32_t>(minimizer_begin.size())};

        for (uint16_t i = 0; i < errors; ++i)
        {
            if (minimizer_begin.size() > 0)
            {
                compute_coverage();
                auto max = std::max_element(coverage.begin(), coverage.end());
                if (i == errors - 1)
                    destroyed += *max;
                else
                {
                    destroyed += *max;
                    std::vector<uint64_t> newBegin;
                    std::vector<uint64_t> newEnd;
                    newBegin.reserve(available - destroyed);
                    newEnd.reserve(available - destroyed);

                    auto idx = std::distance(coverage.begin(), max);
                    auto cb = coverage_begin[idx];
                    auto ce = coverage_end[idx];
                    for (uint64_t i = 0; i < minimizer_begin.size(); ++i)
                    {
                        auto mb = minimizer_begin[i];
                        auto me = minimizer_end[i];
                        if ((mb >= cb && mb <= ce) || (me >= cb && me <= ce))
                            continue;
                        newBegin.push_back(mb);
                        newEnd.push_back(me);
                    }
                    minimizer_begin = std::move(newBegin);
                    minimizer_end = std::move(newEnd);
                    coverage_begin.clear();
                    coverage_end.clear();
                    coverage.clear();
                }
            }
        }
        return destroyed > available ? 0 : available - destroyed;
    }
};

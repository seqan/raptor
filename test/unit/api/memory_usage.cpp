// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h> // for Test, Message, AssertionResult, TestInfo, TestPartResult

#include <string> // for basic_string, string

#include <raptor/argument_parsing/memory_usage.hpp> // for formatted_peak_ram, peak_ram_in_KiB

#if __has_include(<sys/resource.h>)
TEST(peak_ram, as_KiB)
{
    long const result = raptor::peak_ram_in_KiB();
    EXPECT_GT(result, 1024L);
    EXPECT_LT(result, 102400L);
}

TEST(peak_ram, as_string)
{
    std::string const result = raptor::formatted_peak_ram();
    EXPECT_TRUE(result.starts_with("[KiB]: ") || result.starts_with("[MiB]: ")) << result;
}
#endif

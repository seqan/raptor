// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <raptor/argument_parsing/memory_usage.hpp>

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

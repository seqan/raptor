// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <raptor/argument_parsing/memory_usage.hpp>

consteval size_t exact(int const shift)
{
    if (shift >= 64)
        throw std::invalid_argument("Shift value is too big. This is UB!");
    return size_t{1u} << shift;
}

consteval size_t with_decimal(int const shift)
{
    if (shift < 1)
        throw std::invalid_argument("Shift value must be greater than 1.");
    if (shift >= 64)
        throw std::invalid_argument("Shift value is too big. This is UB!");
    return (size_t{1u} << shift) + (size_t{1u} << (shift - 1));
}

consteval size_t rounded(int const shift)
{
    if (shift >= 54)
        throw std::invalid_argument("Shift value is too big. This is UB!");
    return (size_t{1u} << (shift + 10)) - 10;
}

TEST(peak_ram, bytes)
{
    EXPECT_EQ(raptor::detail::formatted_peak_ram(128), "[Bytes]: 128");
}

TEST(peak_ram, KiB)
{
    EXPECT_EQ(raptor::detail::formatted_peak_ram(exact(10)), "[KiB]: 1.0");
    EXPECT_EQ(raptor::detail::formatted_peak_ram(with_decimal(10)), "[KiB]: 1.5");
    EXPECT_EQ(raptor::detail::formatted_peak_ram(rounded(10)), "[MiB]: 1.0");
}

TEST(peak_ram, MiB)
{
    EXPECT_EQ(raptor::detail::formatted_peak_ram(exact(20)), "[MiB]: 1.0");
    EXPECT_EQ(raptor::detail::formatted_peak_ram(with_decimal(20)), "[MiB]: 1.5");
    EXPECT_EQ(raptor::detail::formatted_peak_ram(rounded(20)), "[GiB]: 1.0");
}

TEST(peak_ram, GiB)
{
    EXPECT_EQ(raptor::detail::formatted_peak_ram(exact(30)), "[GiB]: 1.0");
    EXPECT_EQ(raptor::detail::formatted_peak_ram(with_decimal(30)), "[GiB]: 1.5");
    EXPECT_EQ(raptor::detail::formatted_peak_ram(rounded(30)), "[TiB]: 1.0");
}

TEST(peak_ram, TiB)
{
    EXPECT_EQ(raptor::detail::formatted_peak_ram(exact(40)), "[TiB]: 1.0");
    EXPECT_EQ(raptor::detail::formatted_peak_ram(with_decimal(40)), "[TiB]: 1.5");
    EXPECT_EQ(raptor::detail::formatted_peak_ram(rounded(40)), "[PiB]: 1.0");
}

TEST(peak_ram, PiB)
{
    EXPECT_EQ(raptor::detail::formatted_peak_ram(exact(50)), "[PiB]: 1.0");
    EXPECT_EQ(raptor::detail::formatted_peak_ram(with_decimal(50)), "[PiB]: 1.5");
    EXPECT_EQ(raptor::detail::formatted_peak_ram(rounded(50)), "[EiB]: 1.0");
}

TEST(peak_ram, EiB)
{
    EXPECT_EQ(raptor::detail::formatted_peak_ram(exact(60)), "[EiB]: 1.0");
    EXPECT_EQ(raptor::detail::formatted_peak_ram(with_decimal(60)), "[EiB]: 1.5");
}

#if __has_include(<sys/resource.h>)
TEST(peak_ram, via_OS)
{
    std::string const result = raptor::formatted_peak_ram();
    EXPECT_TRUE(result.starts_with("[KiB]: ") || result.starts_with("[MiB]: ")) << result;
}
#endif

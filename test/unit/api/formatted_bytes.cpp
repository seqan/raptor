// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <raptor/argument_parsing/formatted_bytes.hpp>

[[nodiscard]] consteval size_t exact(int const shift)
{
    if (shift >= 64)
        throw std::invalid_argument("Shift value is too big. This is UB!");
    return size_t{1u} << shift;
}

[[nodiscard]] consteval size_t with_decimal(int const shift)
{
    if (shift < 1)
        throw std::invalid_argument("Shift value must be greater than 1.");
    if (shift >= 64)
        throw std::invalid_argument("Shift value is too big. This is UB!");
    return (size_t{1u} << shift) + (size_t{1u} << (shift - 1));
}

[[nodiscard]] consteval size_t rounded(int const shift)
{
    if (shift >= 54)
        throw std::invalid_argument("Shift value is too big. This is UB!");
    return (size_t{1u} << (shift + 10)) - 10;
}

TEST(formatted, bytes)
{
    EXPECT_EQ(raptor::formatted_bytes(128), "[Bytes]: 128");
}

TEST(formatted, KiB)
{
    EXPECT_EQ(raptor::formatted_bytes(exact(10)), "[KiB]: 1.0");
    EXPECT_EQ(raptor::formatted_bytes(with_decimal(10)), "[KiB]: 1.5");
    EXPECT_EQ(raptor::formatted_bytes(rounded(10)), "[MiB]: 1.0");
}

TEST(formatted, MiB)
{
    EXPECT_EQ(raptor::formatted_bytes(exact(20)), "[MiB]: 1.0");
    EXPECT_EQ(raptor::formatted_bytes(with_decimal(20)), "[MiB]: 1.5");
    EXPECT_EQ(raptor::formatted_bytes(rounded(20)), "[GiB]: 1.0");
}

TEST(formatted, GiB)
{
    EXPECT_EQ(raptor::formatted_bytes(exact(30)), "[GiB]: 1.0");
    EXPECT_EQ(raptor::formatted_bytes(with_decimal(30)), "[GiB]: 1.5");
    EXPECT_EQ(raptor::formatted_bytes(rounded(30)), "[TiB]: 1.0");
}

TEST(formatted, TiB)
{
    EXPECT_EQ(raptor::formatted_bytes(exact(40)), "[TiB]: 1.0");
    EXPECT_EQ(raptor::formatted_bytes(with_decimal(40)), "[TiB]: 1.5");
    EXPECT_EQ(raptor::formatted_bytes(rounded(40)), "[PiB]: 1.0");
}

TEST(formatted, PiB)
{
    EXPECT_EQ(raptor::formatted_bytes(exact(50)), "[PiB]: 1.0");
    EXPECT_EQ(raptor::formatted_bytes(with_decimal(50)), "[PiB]: 1.5");
    EXPECT_EQ(raptor::formatted_bytes(rounded(50)), "[EiB]: 1.0");
}

TEST(formatted, EiB)
{
    EXPECT_EQ(raptor::formatted_bytes(exact(60)), "[EiB]: 1.0");
    EXPECT_EQ(raptor::formatted_bytes(with_decimal(60)), "[EiB]: 1.5");
}

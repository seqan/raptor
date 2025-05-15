// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/test/expect_throw_msg.hpp>

#include <raptor/argument_parsing/to_bytes.hpp>

enum class unit
{
    none,
    K,
    Ki,
    M,
    Mi,
    G,
    Gi,
    T,
    Ti,
    P,
    Pi,
    E,
    Ei
};

class to_bytes_test : public testing::TestWithParam<unit>
{};

struct param
{
    std::string str{};
    size_t value{};
};

template <typename value_t>
[[nodiscard]] param get_param(value_t value, unit unit)
{
    std::array<char, 20> buffer{};
    auto [ptr, ec] = std::to_chars(buffer.data(), buffer.data() + buffer.size(), value);
    if (ec != std::errc{})
        throw std::invalid_argument{"Invalid number in \"" + std::string{buffer.data()} + "\"."};

    std::string str{buffer.data(), ptr};

    switch (unit)
    {
    case unit::none:
        return {std::move(str), static_cast<size_t>(value)};
    case unit::K:
        str += "K";
        return {std::move(str), static_cast<size_t>(value * 1000)};
    case unit::Ki:
        str += "Ki";
        return {std::move(str), static_cast<size_t>(value * 1024)};
    case unit::M:
        str += "M";
        return {std::move(str), static_cast<size_t>(value * 1000 * 1000)};
    case unit::Mi:
        str += "Mi";
        return {std::move(str), static_cast<size_t>(value * 1024 * 1024)};
    case unit::G:
        str += "G";
        return {std::move(str), static_cast<size_t>(value * 1000 * 1000 * 1000)};
    case unit::Gi:
        str += "Gi";
        return {std::move(str), static_cast<size_t>(value * 1024 * 1024 * 1024)};
    case unit::T:
        str += "T";
        return {std::move(str), static_cast<size_t>(value * 1000 * 1000 * 1000 * 1000)};
    case unit::Ti:
        str += "Ti";
        return {std::move(str), static_cast<size_t>(value * 1024 * 1024 * 1024 * 1024)};
    case unit::P:
        str += "P";
        return {std::move(str), static_cast<size_t>(value * 1000 * 1000 * 1000 * 1000 * 1000)};
    case unit::Pi:
        str += "Pi";
        return {std::move(str), static_cast<size_t>(value * 1024 * 1024 * 1024 * 1024 * 1024)};
    case unit::E:
        str += "E";
        return {std::move(str), static_cast<size_t>(value * 1000 * 1000 * 1000 * 1000 * 1000 * 1000)};
    case unit::Ei:
        str += "Ei";
        return {std::move(str), static_cast<size_t>(value * 1024 * 1024 * 1024 * 1024 * 1024 * 1024)};
    }

    throw std::invalid_argument{"Unknown unit."};
}

TEST_P(to_bytes_test, unit)
{
    auto const unit = GetParam();
    param param{};

    param = get_param(0ULL, unit);
    EXPECT_EQ(raptor::to_bytes(param.str), param.value);
    param = get_param(0.0, unit);
    EXPECT_EQ(raptor::to_bytes(param.str), param.value);

    param = get_param(1ULL, unit);
    EXPECT_EQ(raptor::to_bytes(param.str), param.value);
    EXPECT_EQ(raptor::to_bytes(" " + param.str + " "), param.value);
    param = get_param(1.0, unit);
    EXPECT_EQ(raptor::to_bytes(param.str), param.value);
    EXPECT_EQ(raptor::to_bytes("  " + param.str + "  "), param.value);

    param = get_param(0.5, unit);
    if (unit == unit::none)
    {
        EXPECT_THROW_MSG((void)raptor::to_bytes(param.str), std::invalid_argument, "Fractional bytes are not allowed.");
    }
    else
    {
        EXPECT_EQ(raptor::to_bytes(param.str), param.value);
    }

    param = get_param(1.5, unit);
    if (unit == unit::none)
    {
        EXPECT_THROW_MSG((void)raptor::to_bytes(param.str), std::invalid_argument, "Fractional bytes are not allowed.");
    }
    else
    {
        EXPECT_EQ(raptor::to_bytes(param.str), param.value);
    }

    param = get_param(1ULL, unit);
    EXPECT_THROW_MSG((void)raptor::to_bytes("-" + param.str),
                     std::invalid_argument,
                     "Negative numbers are not allowed.");
    param = get_param(-1.0, unit);
    EXPECT_THROW_MSG((void)raptor::to_bytes(param.str), std::invalid_argument, "Negative numbers are not allowed.");

    param = get_param(1, unit);
    EXPECT_THROW_MSG((void)raptor::to_bytes(param.str + 'b'),
                     std::invalid_argument,
                     "Unknown unit in \"" + param.str + "b\".");
    EXPECT_THROW_MSG((void)raptor::to_bytes(param.str + 'B'),
                     std::invalid_argument,
                     "Unknown unit in \"" + param.str + "B\".");
}

INSTANTIATE_TEST_SUITE_P(general,
                         to_bytes_test,
                         testing::Values(unit::none,
                                         unit::K,
                                         unit::Ki,
                                         unit::M,
                                         unit::Mi,
                                         unit::G,
                                         unit::Gi,
                                         unit::T,
                                         unit::Ti,
                                         unit::P,
                                         unit::Pi,
                                         unit::E,
                                         unit::Ei),
                         [](testing::TestParamInfo<to_bytes_test::ParamType> const & info)
                         {
                             switch (info.param)
                             {
                             case unit::none:
                                 return "none";
                             case unit::K:
                                 return "K";
                             case unit::Ki:
                                 return "Ki";
                             case unit::M:
                                 return "M";
                             case unit::Mi:
                                 return "Mi";
                             case unit::G:
                                 return "G";
                             case unit::Gi:
                                 return "Gi";
                             case unit::T:
                                 return "T";
                             case unit::Ti:
                                 return "Ti";
                             case unit::P:
                                 return "P";
                             case unit::Pi:
                                 return "Pi";
                             case unit::E:
                                 return "E";
                             case unit::Ei:
                                 return "Ei";
                             }

                             throw std::invalid_argument{"Unknown unit."};
                         });

TEST(to_bytes, empty)
{
    EXPECT_EQ(raptor::to_bytes(""), 0);
    EXPECT_EQ(raptor::to_bytes(" "), 0);
    EXPECT_EQ(raptor::to_bytes("  "), 0);
    EXPECT_EQ(raptor::to_bytes("   "), 0);
}

TEST(to_bytes, no_space_before_unit)
{
    EXPECT_THROW_MSG((void)raptor::to_bytes("100 T"), std::invalid_argument, "Unknown unit in \"100 T\".");
}

TEST(to_bytes, invalid)
{
    EXPECT_THROW_MSG((void)raptor::to_bytes("foo"), std::invalid_argument, "Invalid number in \"foo\".");
}

TEST(to_bytes, out_of_range)
{
    std::string_view str{"1234567890123456789012"};
    EXPECT_THROW_MSG((void)raptor::to_bytes(str),
                     std::out_of_range,
                     "Number in \"" + std::string{str} + "\" is out of range.");
}

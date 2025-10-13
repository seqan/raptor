// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h> // for Test, Message, CmpHelperEQ, CmpHelperEQFailure

#include <cstddef>    // for size_t
#include <filesystem> // for path
#include <string>     // for basic_string, string

#include <raptor/argument_parsing/formatted_index_size.hpp> // for formatted_index_size, index_size_in_KiB

struct helper
{
    static std::filesystem::path data(std::string const & filename)
    {
        return std::filesystem::path{std::string{DATADIR}}.concat(filename);
    }
};

TEST(index_size, ibf)
{
    size_t const size_in_KiB = raptor::index_size_in_KiB(helper::data("1bins23window.index"), 1u);
    EXPECT_EQ(size_in_KiB, 7u);

    std::string const formatted_size = raptor::formatted_index_size(helper::data("1bins23window.index"), 1u);
    EXPECT_EQ(formatted_size, "[KiB]: 7.0");
}

TEST(index_size, ibf_partitioned)
{
    size_t const size_in_KiB = raptor::index_size_in_KiB(helper::data("3.0.partitioned.index"), 4u);
    EXPECT_EQ(size_in_KiB, 56u);

    std::string const formatted_size = raptor::formatted_index_size(helper::data("3.0.partitioned.index"), 4u);
    EXPECT_EQ(formatted_size, "[KiB]: 56.0");
}

TEST(index_size, hibf)
{
    size_t const size_in_KiB = raptor::index_size_in_KiB(helper::data("128bins23window.hibf"), 1u);
    EXPECT_EQ(size_in_KiB, 243u);

    std::string const formatted_size = raptor::formatted_index_size(helper::data("128bins23window.hibf"), 1u);
    EXPECT_EQ(formatted_size, "[KiB]: 243.0");
}

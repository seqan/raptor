// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <raptor/argument_parsing/formatted_index_size.hpp>

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
    size_t const size_in_KiB = raptor::index_size_in_KiB(helper::data("2.0.partitioned.index"), 4u);
    EXPECT_EQ(size_in_KiB, 34u);

    std::string const formatted_size = raptor::formatted_index_size(helper::data("2.0.partitioned.index"), 4u);
    EXPECT_EQ(formatted_size, "[KiB]: 34.0");
}

TEST(index_size, hibf)
{
    size_t const size_in_KiB = raptor::index_size_in_KiB(helper::data("128bins23window.hibf"), 1u);
    EXPECT_EQ(size_in_KiB, 218u);

    std::string const formatted_size = raptor::formatted_index_size(helper::data("128bins23window.hibf"), 1u);
    EXPECT_EQ(formatted_size, "[KiB]: 218.0");
}

// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <raptor/argument_parsing/compute_bin_size.hpp>
#include <raptor/test/cli_test.hpp>

struct compute_bin_size : public raptor_base
{};

TEST_F(compute_bin_size, single_file)
{
    raptor::build_arguments const config{.bin_path = {{data("bin1.fa")}}};
    size_t const result = raptor::compute_bin_size(config);
    EXPECT_EQ(result, 3011u); // exact count = 381, exact size = 3011
}

TEST_F(compute_bin_size, multiple_files)
{
    raptor::build_arguments const config{.bin_path = {{data("bin1.fa"), data("bin2.fa")}}};
    size_t const result = raptor::compute_bin_size(config);
    EXPECT_EQ(result, 3794u); // exact count = 480, exact size = 3794
}

TEST_F(compute_bin_size, multiple_records)
{
    raptor::build_arguments const config{.bin_path = {{data("multi_record_bin.fa")}}};
    size_t const result = raptor::compute_bin_size(config);
    EXPECT_EQ(result, 3794u); // exact count = 480, exact size = 3794
}

TEST_F(compute_bin_size, multiple_files_multiple_records)
{
    raptor::build_arguments const config{.bin_path = {{data("multi_record_bin.fa"), data("multi_record_bin.fa")}}};
    size_t const result = raptor::compute_bin_size(config);
    EXPECT_EQ(result, 3794u); // exact count = 480, exact size = 3794
}

TEST_F(compute_bin_size, mixed)
{
    raptor::build_arguments const config{
        .bin_path = {{data("multi_record_bin.fa")}, {data("bin3.fa")}, {data("bin4.fa")}}};
    size_t const result = raptor::compute_bin_size(config);
    EXPECT_EQ(result, 3794u); // exact count = 480, exact size = 3794
}

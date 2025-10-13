// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h> // for Test, Message, CmpHelperNE, TestPartResult, CmpHelperEQ

#include <algorithm> // for find_if
#include <array>     // for array
#include <charconv>  // for from_chars
#include <ostream>   // for operator<<
#include <string>    // for operator+, basic_string, to_string, operator==, string

#include <sharg/config.hpp>           // for config
#include <sharg/detail/to_string.hpp> // for to_string
#include <sharg/parser.hpp>           // for parser

#include <raptor/argument_parsing/build_arguments.hpp> // for build_arguments
#include <raptor/argument_parsing/shared.hpp>          // for validate_shape

TEST(validate_shape, set_window_to_kmer)
{
    raptor::build_arguments arguments{};
    std::string argv0{"./Test"}, argv1{"--kmer"}, argv2{"30"};
    std::array argv{argv0.data(), argv1.data(), argv2.data()};
    sharg::parser parser{"Test", argv.size(), argv.data()};
    parser.add_option(arguments.kmer_size, sharg::config{.long_id = "kmer"});
    parser.add_option(arguments.window_size, sharg::config{.long_id = "window"});
    parser.add_option(arguments.kmer_size, sharg::config{.long_id = "shape"});

    EXPECT_NE(arguments.kmer_size, 30);
    EXPECT_NE(arguments.window_size, 30);
    parser.parse();
    EXPECT_EQ(arguments.kmer_size, 30);
    EXPECT_NE(arguments.window_size, 30);
    raptor::validate_shape(parser, arguments);
    EXPECT_EQ(arguments.kmer_size, 30);
    EXPECT_EQ(arguments.window_size, 30);
}

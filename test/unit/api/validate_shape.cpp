// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <raptor/argument_parsing/shared.hpp>

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

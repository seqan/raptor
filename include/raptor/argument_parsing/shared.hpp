// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::validate_shape.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <sharg/parser.hpp>

#include <raptor/argument_parsing/build_arguments.hpp>
#include <raptor/argument_parsing/prepare_arguments.hpp>

namespace raptor
{

template <typename argument_t>
    requires std::same_as<argument_t, build_arguments> || std::same_as<argument_t, prepare_arguments>
void validate_shape(sharg::parser & parser, argument_t & arguments)
{
    if (parser.is_option_set("shape"))
    {
        if (parser.is_option_set("kmer"))
            throw sharg::parser_error{"You cannot set both shape and k-mer arguments."};

        uint64_t tmp{};

        std::from_chars(arguments.shape_string.data(),
                        arguments.shape_string.data() + arguments.shape_string.size(),
                        tmp,
                        2);

        arguments.shape = seqan3::shape{seqan3::bin_literal{tmp}};
    }
    else
    {
        arguments.shape = seqan3::shape{seqan3::ungapped{arguments.kmer_size}};
    }

    if (!parser.is_option_set("window"))
        arguments.window_size = arguments.shape.size();
    else if (arguments.shape.size() > arguments.window_size)
        throw sharg::parser_error{"The k-mer size cannot be bigger than the window size."};
}

} // namespace raptor

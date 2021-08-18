// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <raptor/argument_parsing/shared.hpp>
#include <raptor/build/store_index.hpp>

namespace raptor
{

template <bool compressed>
void upgrade_index(upgrade_arguments const & arguments)
{
    constexpr seqan3::data_layout layout = compressed ? seqan3::data_layout::compressed : seqan3::data_layout::uncompressed;
    seqan3::interleaved_bloom_filter<layout> original_index{};

    if (arguments.parts == 1u)
    {
        std::ifstream is{arguments.in_file, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        iarchive(original_index);
        store_index(arguments.out_file, std::move(original_index), arguments);
    }
    else
    {
        for (size_t part : std::views::iota(0u, arguments.parts))
        {
            std::filesystem::path in_file{arguments.in_file};
            in_file += "_" + std::to_string(part);

            std::ifstream is{in_file, std::ios::binary};
            cereal::BinaryInputArchive iarchive{is};
            iarchive(original_index);

            std::filesystem::path out_file{arguments.out_file};
            out_file += "_" + std::to_string(part);
            store_index(out_file, std::move(original_index), arguments);
        }
    }
}

} // namespace raptor

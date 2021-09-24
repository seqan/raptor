// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <string>

namespace raptor::hibf
{

struct build_config
{
    std::string chopper_split_filename{};
    std::string chopper_pack_filename{};
    std::string binning_filename{};
    std::string output_filename{"out.chopper"};
    uint8_t threads{1};
    bool verbose{false};

    size_t hash_funs{2};
    double FPR{0.0001};

    uint8_t k{25};
    uint16_t overlap{200}; // overlap when inserting sequence regions into the IBF
};

} // namespace raptor::hibf

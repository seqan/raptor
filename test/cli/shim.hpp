// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#ifndef __cpp_lib_ranges

#include <range/v3/algorithm/replace.hpp>

namespace std::ranges
{

namespace
{

using ::ranges::cpp20::replace;

} // anonymous namespace

} // namespace std::ranges

#endif // __cpp_lib_ranges

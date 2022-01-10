// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <raptor/search/detail/pascal_row.hpp>

namespace raptor::detail
{

std::vector<size_t> pascal_row(size_t const n)
{
    std::vector<size_t> result(n + 1);
    result[0] = 1;

    for (size_t i = 1; i <= n; ++i)
        result[i] = result[i - 1] * (n + 1 - i) / i;

    return result;
}

} // namespace raptor::detail

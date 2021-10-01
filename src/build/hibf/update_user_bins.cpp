// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <lemon/list_graph.h> // Must be first include.

#include <seqan3/utility/views/join_with.hpp>
#include <seqan3/utility/views/to.hpp>

#include <raptor/build/hibf/update_user_bins.hpp>

namespace raptor::hibf
{

void update_user_bins(build_data & data, std::vector<int64_t> & filename_indices, chopper_pack_record const & record)
{
    size_t const idx = data.request_user_bin_idx();
    data.hibf.user_bins.filename_of_user_bin(idx) = record.filenames | seqan3::views::join_with(std::string{";"}) | seqan3::views::to<std::string>;
    std::fill_n(filename_indices.begin() + record.bin_indices.back(), record.number_of_bins.back(), idx);
}

} // namespace raptor::hibf

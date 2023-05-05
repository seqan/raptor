// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Implements raptor::hibf::update_user_bins.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <lemon/list_graph.h> // Must be first include.

#include <raptor/build/hibf/update_user_bins.hpp>

namespace raptor::hibf
{

void update_user_bins(build_data & data, std::vector<int64_t> & filename_indices, chopper_pack_record const & record)
{
    size_t const idx = data.request_user_bin_idx();

    std::string & user_bin_filenames = data.hibf.user_bins.filename_of_user_bin(idx);
    for (auto const & filename : record.filenames)
    {
        user_bin_filenames += filename;
        user_bin_filenames += ';';
    }
    assert(!user_bin_filenames.empty());
    user_bin_filenames.pop_back();

    std::fill_n(filename_indices.begin() + record.user_bin_info.storage_TB_id,
                record.user_bin_info.number_of_technical_bins,
                idx);
}

} // namespace raptor::hibf

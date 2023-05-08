// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Implements raptor::hibf::chopper_build.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <lemon/list_graph.h> /// Must be first include.

#include <raptor/build/hibf/chopper_build.hpp>
#include <raptor/build/hibf/create_ibfs_from_chopper_pack.hpp>
#include <raptor/build/store_index.hpp>
#include <raptor/index.hpp>

namespace raptor::hibf
{

void chopper_build(build_arguments const & arguments)
{
    build_data data{};

    create_ibfs_from_chopper_pack(data, arguments);

    arguments.index_allocation_timer.start();
    raptor_index<hierarchical_interleaved_bloom_filter> index{window{arguments.window_size},
                                                              arguments.shape,
                                                              arguments.parts,
                                                              data.filenames,
                                                              arguments.fpr,
                                                              std::move(data.hibf)};
    arguments.index_allocation_timer.stop();

    arguments.store_index_timer.start();
    store_index(arguments.out_path, std::move(index), arguments);
    arguments.store_index_timer.stop();
}

} // namespace raptor::hibf

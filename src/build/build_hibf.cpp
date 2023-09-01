// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Implements raptor::build_hibf.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <raptor/build/build_hibf.hpp>
#include <raptor/build/store_index.hpp>
#include <raptor/file_reader.hpp>

#include <hibf/hierarchical_interleaved_bloom_filter.hpp>

namespace raptor
{

void build_hibf(build_arguments const & arguments)
{
    std::variant<file_reader<file_types::sequence>, file_reader<file_types::minimiser>> reader;
    if (arguments.input_is_minimiser)
        reader = file_reader<file_types::minimiser>{};
    else
        reader = file_reader<file_types::sequence>{arguments.shape, arguments.window_size};

    auto input_lambda = [&arguments, &reader](size_t const user_bin_id, seqan::hibf::insert_iterator it)
    {
        std::visit(
            [&](auto const & reader)
            {
                reader.hash_into(arguments.bin_path[user_bin_id], it);
            },
            reader);
    };

    std::ifstream layout_stream{arguments.bin_file};

    seqan::hibf::hierarchical_interleaved_bloom_filter hibf{input_lambda, layout_stream, arguments.threads};

    arguments.index_allocation_timer = std::move(hibf.index_allocation_timer);
    arguments.user_bin_io_timer = std::move(hibf.user_bin_io_timer);
    arguments.merge_kmers_timer = std::move(hibf.merge_kmers_timer);
    arguments.fill_ibf_timer = std::move(hibf.fill_ibf_timer);

    arguments.index_allocation_timer.start();
    raptor_index<index_structure::hibf> index{window{arguments.window_size},
                                              arguments.shape,
                                              arguments.parts,
                                              arguments.bin_path,
                                              arguments.fpr,
                                              std::move(hibf)};
    arguments.index_allocation_timer.stop();

    arguments.store_index_timer.start();
    store_index(arguments.out_path, std::move(index));
    arguments.store_index_timer.stop();
}

} // namespace raptor

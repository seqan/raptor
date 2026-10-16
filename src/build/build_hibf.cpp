// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::build_hibf.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <hibf/hierarchical_interleaved_bloom_filter.hpp>

#include <raptor/build/build_hibf.hpp>
#include <raptor/build/store_index.hpp>
#include <raptor/file_reader.hpp>

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
    std::string const output_path{arguments.out_path};

    assert(arguments.parts >= 1);
    for (size_t part{0}; part < arguments.parts; ++part)
    {
        if (!layout_stream.good())
            throw std::runtime_error{"Expected " + std::to_string(arguments.parts) +
                                     " layouts (given config) but there are only " + std::to_string{part + 1} +
                                     " layouts in file " + arguments.bin_file.string()};

        auto [bin_paths, config, layout] = chopper::read_layout_file(layout_stream);

        // update config and build_arguments
        config.threads = arguments.threads;
        config.input_fn = input_lambda;

        // Inc ase of partitioned HIBF, replace out_path by appending `_[part]` for each index partition.
        arguments.out_path = (arguments.parts > 1) ? output_path + "_" + std::to_string(part) : output_path;

        // construct HIBF
        seqan::hibf::hierarchical_interleaved_bloom_filter hibf{config, layout};

        arguments.index_allocation_timer = std::move(hibf.index_allocation_timer);
        arguments.user_bin_io_timer = std::move(hibf.user_bin_io_timer);
        arguments.merge_kmers_timer = std::move(hibf.merge_kmers_timer);
        arguments.fill_ibf_timer = std::move(hibf.fill_ibf_timer);

        arguments.index_allocation_timer.start();
        raptor_index<index_structure::hibf> index{window{arguments.window_size},
                                                  arguments.shape,
                                                  arguments.parts,
                                                  arguments.bin_path,
                                                  config,
                                                  std::move(hibf)};
        arguments.index_allocation_timer.stop();

        arguments.store_index_timer.start();
        store_index(arguments.out_path, std::move(index));
        arguments.store_index_timer.stop();
    }
}

} // namespace raptor

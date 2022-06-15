// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <seqan3/search/views/minimiser_hash.hpp>

#include <raptor/adjust_seed.hpp>
#include <raptor/build/call_parallel_on_bins.hpp>
#include <raptor/dna4_traits.hpp>
#include <raptor/index.hpp>

namespace raptor
{

template <bool compressed>
class index_factory
{
public:
    index_factory() = default;
    index_factory(index_factory const &) = default;
    index_factory(index_factory &&) = default;
    index_factory & operator=(index_factory const &) = default;
    index_factory & operator=(index_factory &&) = default;
    ~index_factory() = default;

    explicit index_factory(build_arguments const & args) : arguments{std::addressof(args)}
    {}

    template <typename view_t = int>
    [[nodiscard]] auto operator()(view_t && hash_filter_view = 0) const
    {
        auto tmp = construct(std::move(hash_filter_view));

        if constexpr (!compressed)
            return tmp;
        else
            return raptor_index<index_structure::ibf_compressed>{std::move(tmp)};
    }

private:
    build_arguments const * const arguments{nullptr};

    template <typename view_t>
    auto construct(view_t && hash_filter_view) const
    {
        using sequence_file_t = seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>>;

        assert(arguments != nullptr);

        raptor_index<> index{*arguments};

        auto hash_view = [&]()
        {
            if constexpr (std::same_as<view_t, int>)
            {
                return seqan3::views::minimiser_hash(arguments->shape,
                                                     seqan3::window_size{arguments->window_size},
                                                     seqan3::seed{adjust_seed(arguments->shape.count())});
            }
            else
            {
                return seqan3::views::minimiser_hash(arguments->shape,
                                                     seqan3::window_size{arguments->window_size},
                                                     seqan3::seed{adjust_seed(arguments->shape.count())})
                     | hash_filter_view;
            }
        };

        auto worker = [&](auto && zipped_view, auto &&)
        {
            auto & ibf = index.ibf();

            for (auto && [file_names, bin_number] : zipped_view)
                for (auto && file_name : file_names)
                    for (auto && [seq] : sequence_file_t{file_name})
                        for (auto && value : seq | hash_view())
                            ibf.emplace(value, seqan3::bin_index{bin_number});
        };

        call_parallel_on_bins(worker, *arguments);

        return index;
    }
};

} // namespace raptor

// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <seqan3/search/views/minimiser_hash.hpp>

#include <raptor/adjust_seed.hpp>
#include <raptor/build/partition_config.hpp>
#include <raptor/call_parallel_on_bins.hpp>
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

    explicit index_factory(build_arguments const & args, partition_config const & cfg) :
        arguments{std::addressof(args)},
        config{std::addressof(cfg)}
    {}

    [[nodiscard]] auto operator()(size_t const part = 0u) const
    {
        auto tmp = construct(part);

        if constexpr (!compressed)
            return tmp;
        else
            return raptor_index<index_structure::ibf_compressed>{std::move(tmp)};
    }

private:
    build_arguments const * const arguments{nullptr};
    partition_config const * const config{nullptr};

    auto construct(size_t const part) const
    {
        assert(arguments != nullptr);

        raptor_index<> index{*arguments};

        auto worker = [&](auto && zipped_view, auto &&)
        {
            uint64_t hash;
            auto & ibf = index.ibf();

            for (auto && [file_names, bin_number] : zipped_view)
            {
                for (auto && file_name : file_names)
                {
                    std::ifstream infile{file_name, std::ios::binary};
                    if (config == nullptr)
                    {
                        while (infile.read(reinterpret_cast<char *>(&hash), sizeof(hash)))
                            ibf.emplace(hash, seqan3::bin_index{bin_number});
                    }
                    else
                    {
                        while (infile.read(reinterpret_cast<char *>(&hash), sizeof(hash)))
                            if ((hash & config->mask) / config->suffixes_per_part == part)
                                ibf.emplace(hash, seqan3::bin_index{bin_number});
                    }
                }
            }
        };

        call_parallel_on_bins(worker, arguments->bin_path, arguments->threads);

        return index;
    }
};

} // namespace raptor

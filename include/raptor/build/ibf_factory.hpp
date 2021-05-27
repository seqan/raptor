#pragma once

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

#include <raptor/build/call_parallel_on_bins.hpp>

namespace raptor
{

template <bool compressed>
class ibf_factory
{
public:
    ibf_factory(build_arguments const & args) : arguments{args} {}

    ibf_factory() = default;
    ibf_factory(ibf_factory const &) = default;
    ibf_factory(ibf_factory &&) = default;
    ibf_factory & operator=(ibf_factory const &) = default;
    ibf_factory & operator=(ibf_factory &&) = default;
    ~ibf_factory() = default;

    template <typename view_t = std::ranges::empty_view<int>>
        requires (!compressed)
    auto ibf(view_t && hash_filter_view = std::ranges::empty_view<int>()) const
    {
        return construct(std::move(hash_filter_view));
    }

    template <typename view_t = std::ranges::empty_view<int>>
        requires compressed
    auto ibf(view_t && hash_filter_view = std::ranges::empty_view<int>()) const
    {
        auto tmp = construct(std::move(hash_filter_view));

        return seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed>{std::move(tmp)};
    }

private:
    build_arguments const arguments{};

    template <typename view_t = std::ranges::empty_view<int>>
    auto construct(view_t && hash_filter_view = std::ranges::empty_view<int>()) const
    {
        using sequence_file_t = seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>>;

        seqan3::interleaved_bloom_filter<> ibf{seqan3::bin_count{arguments.bins},
                                               seqan3::bin_size{arguments.bits / arguments.parts},
                                               seqan3::hash_function_count{arguments.hash}};

        if constexpr (std::same_as<view_t, std::ranges::empty_view<int>>)
        {
            auto worker = [&] (auto && zipped_view, auto &&)
            {
                auto hash_view = seqan3::views::minimiser_hash(seqan3::ungapped{arguments.kmer_size},
                                                               seqan3::window_size{arguments.window_size},
                                                               seqan3::seed{adjust_seed(arguments.kmer_size)});

                for (auto && [file_name, bin_number] : zipped_view)
                    for (auto && [seq] : sequence_file_t{file_name})
                        for (auto && value : seq | hash_view)
                            ibf.emplace(value, seqan3::bin_index{bin_number});
            };

            call_parallel_on_bins(worker, arguments);
        }
        else
        {
            auto worker = [&] (auto && zipped_view, auto &&)
            {
                auto hash_view = seqan3::views::minimiser_hash(seqan3::ungapped{arguments.kmer_size},
                                                               seqan3::window_size{arguments.window_size},
                                                               seqan3::seed{adjust_seed(arguments.kmer_size)})
                                 | hash_filter_view;

                for (auto && [file_name, bin_number] : zipped_view)
                    for (auto && [seq] : sequence_file_t{file_name})
                        for (auto && value : seq | hash_view)
                            ibf.emplace(value, seqan3::bin_index{bin_number});
            };

            call_parallel_on_bins(worker, arguments);
        }

        return ibf;
    }
};

} // namespace raptor

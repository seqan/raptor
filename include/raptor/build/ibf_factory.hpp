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
    ibf_factory() = default;
    ibf_factory(ibf_factory const &) = default;
    ibf_factory(ibf_factory &&) = default;
    ibf_factory & operator=(ibf_factory const &) = default;
    ibf_factory & operator=(ibf_factory &&) = default;
    ~ibf_factory() = default;

    explicit ibf_factory(build_arguments const & args) : arguments{std::addressof(args)} {}

    template <typename view_t = int>
    [[nodiscard]] auto operator()(view_t && hash_filter_view = 0) const
    {
        auto tmp = construct(std::move(hash_filter_view));

        if constexpr (!compressed)
            return tmp;
        else
            return seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed>{std::move(tmp)};
    }

private:
    build_arguments const * const arguments{nullptr};

    template <typename view_t>
    auto construct(view_t && hash_filter_view) const
    {
        using sequence_file_t = seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>>;

        assert(arguments != nullptr);

        seqan3::interleaved_bloom_filter<> ibf{seqan3::bin_count{arguments->bins},
                                               seqan3::bin_size{arguments->bits / arguments->parts},
                                               seqan3::hash_function_count{arguments->hash}};

        auto hash_view = [&] ()
        {
            if constexpr (!std::ranges::view<view_t>)
            {
                return seqan3::views::minimiser_hash(seqan3::ungapped{arguments->kmer_size},
                                                     seqan3::window_size{arguments->window_size},
                                                     seqan3::seed{adjust_seed(arguments->kmer_size)});
            }
            else
            {
                return seqan3::views::minimiser_hash(seqan3::ungapped{arguments->kmer_size},
                                                     seqan3::window_size{arguments->window_size},
                                                     seqan3::seed{adjust_seed(arguments->kmer_size)})
                       | hash_filter_view;
            }
        };

        auto worker = [&] (auto && zipped_view, auto &&)
        {
            for (auto && [file_names, bin_number] : zipped_view)
                for (auto && file_name : file_names)
                    for (auto && [seq] : sequence_file_t{file_name})
                        for (auto && value : seq | hash_view())
                            ibf.emplace(value, seqan3::bin_index{bin_number});
        };

        call_parallel_on_bins(worker, *arguments);

        return ibf;
    }
};

} // namespace raptor

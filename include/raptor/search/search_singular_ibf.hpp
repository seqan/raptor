// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::search_singular_ibf.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <algorithm>    // for __move, __shuffle, move, shuffle
#include <array>        // for array
#include <cassert>      // for assert
#include <charconv>     // for to_chars
#include <cmath>        // for log2
#include <concepts>     // for same_as
#include <cstddef>      // for size_t
#include <cstdint>      // for uint64_t
#include <filesystem>   // for path
#include <future>       // for launch, async
#include <iterator>     // for back_insert_iterator, back_inserter, operator==
#include <limits>       // for numeric_limits
#include <random>       // for mt19937_64
#include <ranges>       // for operator==, reverse_view, transform_view, common_view
#include <span>         // for span
#include <string>       // for basic_string, string
#include <string_view>  // for basic_string_view, string_view
#include <system_error> // for errc
#include <tuple>        // for get
#include <vector>       // for vector

#include <seqan3/alphabet/nucleotide/dna4.hpp>         // for dna4
#include <seqan3/contrib/std/chunk_view.hpp>           // for chunk, chunk_fn, chunk_view, iter_move, operator==
#include <seqan3/contrib/std/detail/adaptor_base.hpp>  // for operator|
#include <seqan3/core/range/detail/adaptor_base.hpp>   // for operator|
#include <seqan3/io/detail/misc.hpp>                   // for set_format
#include <seqan3/io/record.hpp>                        // for fields, field
#include <seqan3/io/sequence_file/input.hpp>           // for sequence_file_input
#include <seqan3/io/sequence_file/record.hpp>          // for sequence_record
#include <seqan3/search/views/kmer_hash.hpp>           // for operator-, operator==
#include <seqan3/search/views/minimiser.hpp>           // for minimiser_view, operator==
#include <seqan3/search/views/minimiser_hash.hpp>      // for seed, window_size, minimiser_hash, minimiser_hash_fn
#include <seqan3/utility/container/dynamic_bitset.hpp> // for operator==

#include <hibf/misc/timer.hpp> // for serial_timer, concurrent_timer

#include <raptor/adjust_seed.hpp>                       // for adjust_seed
#include <raptor/argument_parsing/search_arguments.hpp> // for search_arguments
#include <raptor/dna4_traits.hpp>                       // for dna4_traits
#include <raptor/index.hpp>                             // for ibf, raptor_index
#include <raptor/search/do_parallel.hpp>                // for do_parallel
#include <raptor/search/load_index.hpp>                 // for load_index
#include <raptor/search/sync_out.hpp>                   // for sync_out
#include <raptor/threshold/threshold.hpp>               // for threshold

namespace raptor
{

template <typename index_t>
void search_singular_ibf(search_arguments const & arguments, index_t && index)
{
    constexpr bool is_ibf = std::same_as<index_t, raptor_index<index_structure::ibf>>;

    auto cereal_future = std::async(std::launch::async,
                                    [&]()
                                    {
                                        load_index(index, arguments);
                                    });

    seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::id, seqan3::field::seq>> fin{
        arguments.query_file};
    using record_type = typename decltype(fin)::record_type;
    std::vector<record_type> records{};

    sync_out synced_out{arguments};

    raptor::threshold::threshold const thresholder{arguments.make_threshold_parameters()};

    auto worker = [&](size_t const start, size_t const extent)
    {
        seqan::hibf::serial_timer local_compute_minimiser_timer{};
        seqan::hibf::serial_timer local_query_ibf_timer{};
        seqan::hibf::serial_timer local_generate_results_timer{};

        auto agent = index.ibf().membership_agent();

        std::string result_string{};
        std::array<char, std::numeric_limits<uint64_t>::digits10 + 1> buffer{};
        std::vector<uint64_t> minimiser;

        auto hash_adaptor = seqan3::views::minimiser_hash(arguments.shape,
                                                          seqan3::window_size{arguments.window_size},
                                                          seqan3::seed{adjust_seed(arguments.shape_weight)});

        for (auto && [id, seq] : std::span{records.data() + start, extent})
        {
            result_string.clear();
            result_string += id;
            result_string += '\t';

            auto minimiser_view = seq | hash_adaptor | std::views::common;
            local_compute_minimiser_timer.start();
            minimiser.assign(minimiser_view.begin(), minimiser_view.end());
            local_compute_minimiser_timer.stop();

            size_t const minimiser_count{minimiser.size()};
            size_t const threshold = thresholder.get(minimiser_count);

            local_query_ibf_timer.start();
            auto & user_bin_ids = agent.membership_for(minimiser, threshold);
            local_query_ibf_timer.stop();
            local_generate_results_timer.start();
            for (auto && user_bin : user_bin_ids)
            {
                auto conv = std::to_chars(buffer.data(), buffer.data() + buffer.size(), user_bin);
                assert(conv.ec == std::errc{});
                std::string_view sv{buffer.data(), conv.ptr};
                result_string += sv;
                result_string += ',';
            }

            if (auto & last_char = result_string.back(); last_char == ',')
                last_char = '\n';
            else
                result_string += '\n';

            synced_out.write(result_string);
            local_generate_results_timer.stop();
        }

        arguments.compute_minimiser_timer += local_compute_minimiser_timer;
        arguments.query_ibf_timer += local_query_ibf_timer;
        arguments.generate_results_timer += local_generate_results_timer;
    };

    auto write_header = [&]()
    {
        if constexpr (is_ibf)
            return synced_out.write_header(arguments, index.ibf().hash_function_count());
        else
            return synced_out.write_header(arguments, index.ibf().ibf_vector[0].hash_function_count());
    };

    for (auto && chunked_records : fin | seqan::stl::views::chunk((1ULL << 20) * 10))
    {
        records.clear();
        arguments.query_file_io_timer.start();
        std::ranges::move(chunked_records, std::back_inserter(records));
        // Very fast, improves parallel processing when chunks of the query belong to the same bin.
        std::ranges::shuffle(records, std::mt19937_64{0u});
        arguments.query_file_io_timer.stop();

        cereal_future.get();
        [[maybe_unused]] static bool header_written = write_header(); // called exactly once

        arguments.parallel_search_timer.start();
        do_parallel(worker, records.size(), arguments.threads);
        arguments.parallel_search_timer.stop();
    }
}

} // namespace raptor

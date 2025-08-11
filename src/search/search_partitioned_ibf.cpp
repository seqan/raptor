// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::search_partitioned_ibf.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <algorithm>    // for __move, __shuffle, move, shuffle
#include <array>        // for array
#include <cassert>      // for assert
#include <charconv>     // for to_chars, to_chars_result
#include <cmath>        // for log2
#include <cstddef>      // for size_t
#include <cstdint>      // for uint16_t, uint64_t
#include <filesystem>   // for path
#include <future>       // for async, launch, future
#include <iterator>     // for back_insert_iterator, operator==, back_inserter
#include <limits>       // for numeric_limits
#include <random>       // for mt19937_64
#include <ranges>       // for transform_view, reverse_view, views, __pipeable, ope...
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

#include <hibf/interleaved_bloom_filter.hpp> // for interleaved_bloom_filter
#include <hibf/misc/counting_vector.hpp>     // for counting_vector
#include <hibf/misc/timer.hpp>               // for serial_timer, concurrent_timer

#include <raptor/adjust_seed.hpp>                       // for adjust_seed
#include <raptor/argument_parsing/search_arguments.hpp> // for search_arguments
#include <raptor/build/partition_config.hpp>            // for partition_config
#include <raptor/dna4_traits.hpp>                       // for dna4_traits
#include <raptor/index.hpp>                             // for raptor_index, ibf
#include <raptor/search/do_parallel.hpp>                // for do_parallel
#include <raptor/search/load_index.hpp>                 // for load_index
#include <raptor/search/search_partitioned_ibf.hpp>     // for search_partitioned_ibf
#include <raptor/search/sync_out.hpp>                   // for sync_out
#include <raptor/threshold/threshold.hpp>               // for threshold

namespace raptor
{

void search_partitioned_ibf(search_arguments const & arguments)
{
    auto index = raptor_index<index_structure::ibf>{};
    partition_config const cfg{arguments.parts};

    seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::id, seqan3::field::seq>> fin{
        arguments.query_file};
    using record_type = typename decltype(fin)::record_type;
    std::vector<record_type> records{};

    sync_out synced_out{arguments};

    auto write_header = [&]()
    {
        return synced_out.write_header(arguments, index.ibf().hash_function_count());
    };

    raptor::threshold::threshold const thresholder{arguments.make_threshold_parameters()};

    for (auto && chunked_records : fin | seqan::stl::views::chunk((1ULL << 20) * 10))
    {
        auto cereal_future = std::async(std::launch::async,
                                        [&]() // GCOVR_EXCL_LINE
                                        {
                                            load_index(index, arguments, 0);
                                        });

        records.clear();
        arguments.query_file_io_timer.start();
        std::ranges::move(chunked_records, std::back_inserter(records));
        // Very fast, improves parallel processing when chunks of the query belong to the same bin.
        std::ranges::shuffle(records, std::mt19937_64{0u});
        arguments.query_file_io_timer.stop();

        cereal_future.get();
        [[maybe_unused]] static bool header_written = write_header(); // called exactly once

        std::vector<seqan::hibf::counting_vector<uint16_t>> counts(
            records.size(),
            seqan::hibf::counting_vector<uint16_t>(index.ibf().bin_count(), 0));

        size_t part{};

        auto count_task = [&](size_t const start, size_t const extent)
        {
            seqan::hibf::serial_timer local_compute_minimiser_timer{};
            seqan::hibf::serial_timer local_query_ibf_timer{};

            auto & ibf = index.ibf();
            auto counter = ibf.template counting_agent<uint16_t>();
            size_t counter_id = start;
            std::vector<uint64_t> minimiser;

            auto hash_view = seqan3::views::minimiser_hash(arguments.shape,
                                                           seqan3::window_size{arguments.window_size},
                                                           seqan3::seed{adjust_seed(arguments.shape_weight)});

            for (auto && [id, seq] : std::span{records.data() + start, extent})
            {
                auto minimiser_view = seq | hash_view | std::views::common;
                local_compute_minimiser_timer.start();
                minimiser.assign(minimiser_view.begin(), minimiser_view.end());
                local_compute_minimiser_timer.stop();

                // GCOVR_EXCL_START
                auto filtered = minimiser
                              | std::views::filter(
                                    [&](auto && hash)
                                    {
                                        return cfg.hash_partition(hash) == part;
                                    });
                // GCOVR_EXCL_STOP

                local_query_ibf_timer.start();
                counts[counter_id++] += counter.bulk_count(filtered);
                local_query_ibf_timer.stop();
            }

            arguments.compute_minimiser_timer += local_compute_minimiser_timer;
            arguments.query_ibf_timer += local_query_ibf_timer;
        };

        arguments.parallel_search_timer.start();
        do_parallel(count_task, records.size(), arguments.threads);
        arguments.parallel_search_timer.stop();
        ++part;

        for (; part < arguments.parts - 1u; ++part)
        {
            load_index(index, arguments, part);
            arguments.parallel_search_timer.start();
            do_parallel(count_task, records.size(), arguments.threads);
            arguments.parallel_search_timer.stop();
        }

        assert(part == arguments.parts - 1u);
        load_index(index, arguments, part);

        auto output_task = [&](size_t const start, size_t const extent)
        {
            seqan::hibf::serial_timer local_compute_minimiser_timer{};
            seqan::hibf::serial_timer local_query_ibf_timer{};
            seqan::hibf::serial_timer local_generate_results_timer{};

            auto & ibf = index.ibf();
            auto counter = ibf.template counting_agent<uint16_t>();
            size_t counter_id = start;
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

                // GCOVR_EXCL_START
                auto filtered = minimiser
                              | std::views::filter(
                                    [&](auto && hash)
                                    {
                                        return cfg.hash_partition(hash) == part;
                                    });
                // GCOVR_EXCL_STOP
                local_query_ibf_timer.start();
                counts[counter_id] += counter.bulk_count(filtered);
                local_query_ibf_timer.stop();

                size_t const minimiser_count{minimiser.size()};
                size_t current_bin{0};

                size_t const threshold = thresholder.get(minimiser_count);
                local_generate_results_timer.start();
                for (auto && count : counts[counter_id++])
                {
                    if (count >= threshold)
                    {
                        auto conv = std::to_chars(buffer.data(), buffer.data() + buffer.size(), current_bin);
                        assert(conv.ec == std::errc{});
                        std::string_view sv{buffer.data(), conv.ptr};
                        result_string += sv;
                        result_string += ',';
                    }
                    ++current_bin;
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

        arguments.parallel_search_timer.start();
        do_parallel(output_task, records.size(), arguments.threads);
        arguments.parallel_search_timer.stop();
    }
}

} // namespace raptor

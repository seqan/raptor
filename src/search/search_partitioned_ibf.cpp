// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Implements raptor::search_partitioned_ibf.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <seqan3/search/views/minimiser_hash.hpp>

#include <raptor/adjust_seed.hpp>
#include <raptor/build/partition_config.hpp>
#include <raptor/dna4_traits.hpp>
#include <raptor/search/do_parallel.hpp>
#include <raptor/search/load_index.hpp>
#include <raptor/search/search_partitioned_ibf.hpp>
#include <raptor/search/sync_out.hpp>
#include <raptor/threshold/threshold.hpp>

namespace raptor
{

template <bool compressed>
void search_partitioned_ibf(search_arguments const & arguments)
{
    using index_structure_t = std::conditional_t<compressed, index_structure::ibf_compressed, index_structure::ibf>;
    auto index = raptor_index<index_structure_t>{};
    partition_config const cfg{arguments.parts};

    seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::id, seqan3::field::seq>> fin{
        arguments.query_file};
    using record_type = typename decltype(fin)::record_type;
    std::vector<record_type> records{};

    auto cereal_worker = [&]()
    {
        load_index(index, arguments, 0);
    };

    sync_out synced_out{arguments};

    raptor::threshold::threshold const thresholder{arguments.make_threshold_parameters()};

    for (auto && chunked_records : fin | seqan3::views::chunk((1ULL << 20) * 10))
    {
        auto cereal_handle = std::async(std::launch::async, cereal_worker);

        records.clear();
        arguments.query_file_io_timer.start();
        std::ranges::move(chunked_records, std::back_inserter(records));
        arguments.query_file_io_timer.stop();

        cereal_handle.wait();

        std::vector<seqan3::counting_vector<uint16_t>> counts(
            records.size(),
            seqan3::counting_vector<uint16_t>(index.ibf().bin_count(), 0));

        size_t part{};

        auto count_task = [&](size_t const start, size_t const end)
        {
            timer<concurrent::no> local_compute_minimiser_timer{};
            timer<concurrent::no> local_query_ibf_timer{};

            auto & ibf = index.ibf();
            auto counter = ibf.template counting_agent<uint16_t>();
            size_t counter_id = start;
            std::vector<uint64_t> minimiser;

            auto hash_view = seqan3::views::minimiser_hash(arguments.shape,
                                                           seqan3::window_size{arguments.window_size},
                                                           seqan3::seed{adjust_seed(arguments.shape_weight)});

            for (auto && [id, seq] : records | seqan3::views::slice(start, end))
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

        do_parallel(count_task, records.size(), arguments.threads);
        ++part;

        for (; part < arguments.parts - 1u; ++part)
        {
            load_index(index, arguments, part);
            do_parallel(count_task, records.size(), arguments.threads);
        }

        assert(part == arguments.parts - 1u);
        load_index(index, arguments, part);

        auto output_task = [&](size_t const start, size_t const end)
        {
            timer<concurrent::no> local_compute_minimiser_timer{};
            timer<concurrent::no> local_query_ibf_timer{};
            timer<concurrent::no> local_generate_results_timer{};

            auto & ibf = index.ibf();
            auto counter = ibf.template counting_agent<uint16_t>();
            size_t counter_id = start;
            std::string result_string{};
            std::vector<uint64_t> minimiser;

            auto hash_adaptor = seqan3::views::minimiser_hash(arguments.shape,
                                                              seqan3::window_size{arguments.window_size},
                                                              seqan3::seed{adjust_seed(arguments.shape_weight)});

            for (auto && [id, seq] : records | seqan3::views::slice(start, end))
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
                        result_string += std::to_string(current_bin);
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

        do_parallel(output_task, records.size(), arguments.threads);
    }
}

template void search_partitioned_ibf<false>(search_arguments const & arguments);

template void search_partitioned_ibf<true>(search_arguments const & arguments);

} // namespace raptor

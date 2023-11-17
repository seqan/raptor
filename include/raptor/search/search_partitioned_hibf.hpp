// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides raptor::search_partitioned_hibf.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <seqan3/search/views/minimiser_hash.hpp>

#include <raptor/adjust_seed.hpp>
#include <raptor/contrib/std/chunk_view.hpp>
#include <raptor/dna4_traits.hpp>
#include <raptor/search/do_parallel.hpp>
#include <raptor/search/load_index.hpp>
#include <raptor/search/sync_out.hpp>
#include <raptor/threshold/threshold.hpp>

namespace raptor
{

template <typename index_t>
void search_partitioned_hibf(search_arguments const & arguments, index_t && index)
{
    seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::id, seqan3::field::seq>> fin{
        arguments.query_file};
    using record_type = typename decltype(fin)::record_type;

    std::vector<record_type> records{};

    sync_out synced_out{arguments};

    std::vector<std::string> results; // cache results since we are searching multiple hibfs

    raptor::threshold::threshold const thresholder{arguments.make_threshold_parameters()};

    // searching with storing all results in results map
    auto worker = [&](size_t const start, size_t const extent, bool const output_results)
    {
        seqan::hibf::serial_timer local_compute_minimiser_timer{};
        seqan::hibf::serial_timer local_query_ibf_timer{};
        seqan::hibf::serial_timer local_generate_results_timer{};

        auto counter = index.ibf().membership_agent();

        std::string result_string{};
        std::vector<uint64_t> minimiser;

        auto hash_adaptor = seqan3::views::minimiser_hash(arguments.shape,
                                                          seqan3::window_size{arguments.window_size},
                                                          seqan3::seed{adjust_seed(arguments.shape_weight)});

        for (size_t pos = start; pos < start + extent; ++pos)
        {
            auto const & seq = records[pos].sequence();
            std::string & result_string = results[pos];

            auto minimiser_view = seq | hash_adaptor | std::views::common;
            local_compute_minimiser_timer.start();
            minimiser.assign(minimiser_view.begin(), minimiser_view.end());
            local_compute_minimiser_timer.stop();

            size_t const minimiser_count{minimiser.size()};
            size_t const threshold = thresholder.get(minimiser_count);

            local_query_ibf_timer.start();
            auto & result = counter.membership_for(minimiser, threshold); // Results contains user bin IDs
            local_query_ibf_timer.stop();
            local_generate_results_timer.start();
            for (auto && count : result)
            {
                result_string += std::to_string(count);
                result_string += ',';
            }

            if (output_results)
            {
                result_string.insert(result_string.begin(), '\t');
                auto const & id = records[pos].id();
                result_string.insert(result_string.begin(), id.begin(), id.end());

                if (auto & last_char = result_string.back(); last_char == ',')
                    last_char = '\n';
                else
                    result_string += '\n';

                synced_out.write(result_string);
                result_string.clear(); // free memory
            }
            local_generate_results_timer.stop();
        }

        arguments.compute_minimiser_timer += local_compute_minimiser_timer;
        arguments.query_ibf_timer += local_query_ibf_timer;
        arguments.generate_results_timer += local_generate_results_timer;
    };

    for (auto && chunked_records : fin | seqan::stl::views::chunk((1ULL << 20) * 10))
    {
        // prefetch the first partition while query IO is done
        auto cereal_future = std::async(std::launch::async,
                                        [&]()
                                        {
                                            load_index(index, arguments, 0);
                                        });

        records.clear();
        arguments.query_file_io_timer.start();
        std::ranges::move(chunked_records, std::back_inserter(records));
        arguments.query_file_io_timer.stop();

        results.resize(records.size());

        cereal_future.get();
        synced_out.write_header(arguments, index.ibf().ibf_vector[0].hash_function_count());

        assert(arguments.parts > 0);
        for (int part = 0; part < arguments.parts - 1; ++part)
        {
            do_parallel(worker, records.size(), arguments.threads, false/*do not write results*/);
            load_index(index, arguments, part + 1);
        }

        do_parallel(worker, records.size(), arguments.threads, true/*write results*/);
    }
}

} // namespace raptor

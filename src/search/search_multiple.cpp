// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <seqan3/search/views/minimiser_hash.hpp>

#include <raptor/adjust_seed.hpp>
#include <raptor/dna4_traits.hpp>
#include <raptor/search/do_parallel.hpp>
#include <raptor/search/load_index.hpp>
#include <raptor/search/search_multiple.hpp>
#include <raptor/search/sync_out.hpp>
#include <raptor/search/threshold.hpp>

namespace raptor
{

template <bool compressed>
void search_multiple(search_arguments const & arguments)
{
    using index_structure_t = std::conditional_t<compressed, index_structure::ibf_compressed, index_structure::ibf>;
    auto index = raptor_index<index_structure_t>{};

    seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::id, seqan3::field::seq>> fin{arguments.query_file};
    using record_type = typename decltype(fin)::record_type;
    std::vector<record_type> records{};

    double index_io_time{0.0};
    double reads_io_time{0.0};
    double compute_time{0.0};

    auto cereal_worker = [&] ()
    {
        load_index(index, arguments, 0, index_io_time);
    };

    sync_out synced_out{arguments.out_file};

    {
        size_t position{};
        std::string line{};
        for (auto const & file_list : arguments.bin_path)
        {
            line.clear();
            line = '#';
            line += std::to_string(position);
            line += '\t';
            for (auto const & filename : file_list)
            {
                line += filename;
                line += ',';
            }
            line.back() = '\n';
            synced_out << line;
            ++position;
        }
        synced_out << "#QUERY_NAME\tUSER_BINS\n";
    }

    threshold const thresholder{arguments};

    for (auto && chunked_records : fin | seqan3::views::chunk((1ULL<<20)*10))
    {
        auto cereal_handle = std::async(std::launch::async, cereal_worker);

        records.clear();
        auto start = std::chrono::high_resolution_clock::now();
        std::ranges::move(chunked_records, std::cpp20::back_inserter(records));
        auto end = std::chrono::high_resolution_clock::now();
        reads_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();

        cereal_handle.wait();

        std::vector<seqan3::counting_vector<uint16_t>> counts(records.size(),
                                                              seqan3::counting_vector<uint16_t>(index.ibf().bin_count(), 0));

        auto count_task = [&](size_t const start, size_t const end)
        {
            auto & ibf = index.ibf();
            auto counter = ibf.template counting_agent<uint16_t>();
            size_t counter_id = start;

            auto hash_view = seqan3::views::minimiser_hash(arguments.shape,
                                                           seqan3::window_size{arguments.window_size},
                                                           seqan3::seed{adjust_seed(arguments.shape_weight)});

            for (auto && [id, seq] : records | seqan3::views::slice(start, end))
            {
                (void) id;
                auto & result = counter.bulk_count(seq | hash_view);
                counts[counter_id++] += result;
            }
        };

        do_parallel(count_task, records.size(), arguments.threads, compute_time);

        for (size_t const part : std::views::iota(1u, static_cast<unsigned int>(arguments.parts - 1)))
        {
            load_index(index, arguments, part, index_io_time);
            do_parallel(count_task, records.size(), arguments.threads, compute_time);
        }

        load_index(index, arguments, arguments.parts - 1, index_io_time);

        auto output_task = [&](size_t const start, size_t const end)
        {
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
                minimiser.assign(minimiser_view.begin(), minimiser_view.end());

                counts[counter_id] += counter.bulk_count(minimiser);
                size_t const minimiser_count{minimiser.size()};
                size_t current_bin{0};

                size_t const threshold = thresholder.get(minimiser_count);
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
            }
        };

        do_parallel(output_task, records.size(), arguments.threads, compute_time);
    }

// LCOV_EXCL_START
    if (arguments.write_time)
    {
        std::filesystem::path file_path{arguments.out_file};
        file_path += ".time";
        std::ofstream file_handle{file_path};
        file_handle << "Index I/O\tReads I/O\tCompute\n";
        file_handle << std::fixed
                    << std::setprecision(2)
                    << index_io_time << '\t'
                    << reads_io_time << '\t'
                    << compute_time;
    }
// LCOV_EXCL_STOP
}

template
void search_multiple<false>(search_arguments const & arguments);

template
void search_multiple<true>(search_arguments const & arguments);

} // namespace raptor

#pragma once

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

#include <raptor/search/compute_simple_model.hpp>
#include <raptor/search/do_parallel.hpp>
#include <raptor/search/load_ibf.hpp>
#include <raptor/search/sync_out.hpp>

namespace raptor
{

template <bool compressed>
void run_program_single(search_arguments const & arguments)
{
    constexpr seqan3::data_layout ibf_data_layout = compressed ? seqan3::data_layout::compressed :
                                                                 seqan3::data_layout::uncompressed;
    auto ibf = seqan3::interleaved_bloom_filter<ibf_data_layout>{};

    double ibf_io_time{0.0};
    double reads_io_time{0.0};
    double compute_time{0.0};

    auto cereal_worker = [&] ()
    {
        load_ibf(ibf, arguments, ibf_io_time);
    };
    auto cereal_handle = std::async(std::launch::async, cereal_worker);

    seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::id, seqan3::field::seq>> fin{arguments.query_file};
    using record_type = typename decltype(fin)::record_type;
    std::vector<record_type> records{};

    sync_out synced_out{arguments.out_file};

    size_t const kmers_per_window = arguments.window_size - arguments.kmer_size + 1;
    size_t const kmers_per_pattern = arguments.pattern_size - arguments.kmer_size + 1;
    size_t const min_number_of_minimisers = kmers_per_window == 1 ? kmers_per_pattern :
                                                std::ceil(kmers_per_pattern / static_cast<double>(kmers_per_window));
    size_t const kmer_lemma = arguments.pattern_size + 1u > (arguments.errors + 1u) * arguments.kmer_size ?
                                arguments.pattern_size + 1u - (arguments.errors + 1u) * arguments.kmer_size :
                                0;
    size_t const max_number_of_minimisers = arguments.pattern_size - arguments.window_size + 1;
    std::vector<size_t> const precomp_thresholds = compute_simple_model(arguments);

    auto worker = [&] (size_t const start, size_t const end)
    {
        auto counter = ibf.template counting_agent<uint16_t>();
        std::string result_string{};
        std::vector<uint64_t> minimiser;

        auto hash_view = seqan3::views::minimiser_hash(seqan3::ungapped{arguments.kmer_size},
                                                       seqan3::window_size{arguments.window_size},
                                                       seqan3::seed{adjust_seed(arguments.kmer_size)});

        for (auto && [id, seq] : records | seqan3::views::slice(start, end))
        {
            minimiser.clear();
            result_string.clear();
            result_string += id;
            result_string += '\t';

            minimiser = seq | hash_view | seqan3::views::to<std::vector<uint64_t>>;
            auto & result = counter.bulk_count(minimiser);
            size_t const minimiser_count{minimiser.size()};
            size_t current_bin{0};

            size_t const threshold = arguments.treshold_was_set ?
                                         static_cast<size_t>(minimiser_count * arguments.threshold) :
                                         kmers_per_window == 1 ? kmer_lemma :
                                         precomp_thresholds[std::min(minimiser_count < min_number_of_minimisers ?
                                                                         0 :
                                                                         minimiser_count - min_number_of_minimisers,
                                                                     max_number_of_minimisers -
                                                                         min_number_of_minimisers)] + 2;

            for (auto && count : result)
            {
                if (count >= threshold)
                {
                    result_string += std::to_string(current_bin);
                    result_string += ',';
                }
                ++current_bin;
            }
            result_string += '\n';
            synced_out.write(result_string);
        }
    };

    for (auto && chunked_records : fin | seqan3::views::chunk((1ULL<<20)*10))
    {
        records.clear();
        auto start = std::chrono::high_resolution_clock::now();
        std::ranges::move(chunked_records, std::cpp20::back_inserter(records));
        auto end = std::chrono::high_resolution_clock::now();
        reads_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();

        cereal_handle.wait();

        do_parallel(worker, records.size(), arguments.threads, compute_time);
    }

    if (arguments.write_time)
    {
        std::filesystem::path file_path{arguments.out_file};
        file_path += ".time";
        std::ofstream file_handle{file_path};
        file_handle << "IBF I/O\tReads I/O\tCompute\n";
        file_handle << std::fixed
                    << std::setprecision(2)
                    << ibf_io_time << '\t'
                    << reads_io_time << '\t'
                    << compute_time;
    }
}

} // namespace raptor

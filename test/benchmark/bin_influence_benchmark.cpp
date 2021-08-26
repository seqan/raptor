// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/utility/views/slice.hpp>

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>

static constexpr size_t const genome_size{5000}; // 4'300'000'000
static constexpr size_t const read_size{100};
static constexpr size_t const read_count{1000}; // 1'000'000
static constexpr size_t const ibf_size{8'388'608/*=1MiB*/}; // 34'359'738'368/*=4GiB*/

static std::vector<seqan3::dna4> const genome{seqan3::test::generate_sequence<seqan3::dna4>(genome_size, 0, 0)};
static std::vector<std::vector<seqan3::dna4>> const reads{[] (auto const & genome) {
    std::vector<std::vector<seqan3::dna4>> result(read_count);
    size_t i{};
    for (auto && read_start : seqan3::test::generate_numeric_sequence<size_t>(read_count, 0, genome_size - read_size + 1, 0))
    {
        auto v = genome | seqan3::views::slice(read_start, read_start + read_size);
        result[i++].assign(v.begin(), v.end());
    }
    return result;
    } (genome)};

static void search_benchmark(benchmark::State & state)
{
    size_t const bin_count = static_cast<size_t>(state.range(0));
    size_t const hash_num{2u};
    size_t const bin_size{ibf_size / bin_count};
    size_t const chunk_size{(genome_size + bin_count - 1) / bin_count};

    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf{seqan3::bin_count{bin_count},
                                                                            seqan3::bin_size{bin_size},
                                                                            seqan3::hash_function_count{hash_num}};

    size_t bin_counter{};
    for (auto && sequence : genome | seqan3::views::chunk(chunk_size))
        for (auto && hash : sequence | seqan3::views::kmer_hash(seqan3::ungapped{19u}))
            ibf.emplace(hash, seqan3::bin_index{bin_counter++});

    auto agent = ibf.counting_agent<uint16_t>();
    for (auto _ : state)
        for (auto && query : reads)
            benchmark::DoNotOptimize(agent.bulk_count(query | seqan3::views::kmer_hash(seqan3::ungapped{19u})));
}

BENCHMARK(search_benchmark)->RangeMultiplier(2)->Range(64, 65536);

BENCHMARK_MAIN();

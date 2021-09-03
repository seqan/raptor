// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/algorithm/detail/execution_handler_parallel.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/utility/views/slice.hpp>
#include <seqan3/utility/views/zip.hpp>

#include <raptor/adjust_seed.hpp>

constexpr size_t operator""_MiB(unsigned long long int number)
{
    return number << 23;
}

constexpr size_t operator""_GiB(unsigned long long int number)
{
    return number << 33;
}

#if 1
static constexpr size_t const genome_size{5000};
static constexpr size_t const read_size{100};
static constexpr size_t const read_count{1000};
static constexpr size_t const ibf_size{1_MiB};
#else
static constexpr size_t const genome_size{4'300'000'000};
static constexpr size_t const read_size{100};
static constexpr size_t const read_count{10'000'000};
static constexpr size_t const ibf_size{4_GiB};
#endif

static constexpr seqan3::window_size const window_size{19u};
static constexpr seqan3::shape const shape{seqan3::ungapped{19u}};
static constexpr size_t const construct_threads{1};
using ibf_t = seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>;

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

static ibf_t construct_ibf(size_t const bin_count, auto && hash_adaptor)
{
    size_t const hash_num{2u};
    size_t const bin_size{ibf_size / bin_count};
    size_t const chunk_size{(genome_size + bin_count - 1) / bin_count};
    size_t const chunked_genome_size{(genome_size + bin_count - 1) / bin_count};
    size_t const workload_size = std::clamp<size_t>(std::bit_ceil(bin_count / construct_threads), 8u, 64u);
    ibf_t ibf{seqan3::bin_count{bin_count}, seqan3::bin_size{bin_size}, seqan3::hash_function_count{hash_num}};

    auto chunked_genomes = genome | seqan3::views::chunk(chunked_genome_size);
    auto workload = seqan3::views::zip(chunked_genomes, std::views::iota(0u)) | seqan3::views::chunk(workload_size);
    auto worker = [&ibf, &hash_adaptor] (auto && payload, auto &&)
    {
        for (auto && [sequence, bin_number] : payload)
            for (auto && hash : sequence | hash_adaptor)
                    ibf.emplace(hash, seqan3::bin_index{bin_number});
    };

    seqan3::detail::execution_handler_parallel executioner{construct_threads};
    executioner.bulk_execute(std::move(worker), std::move(workload), [](){});

    return ibf;
}

static void bulk_count(benchmark::State & state)
{
    size_t const bin_count = static_cast<size_t>(state.range(0));
    auto hash_adaptor = seqan3::views::minimiser_hash(shape,
                                                      window_size,
                                                      seqan3::seed{raptor::adjust_seed(shape.count())});

    ibf_t const ibf{construct_ibf(bin_count, hash_adaptor)};

    auto agent = ibf.counting_agent<uint16_t>();
    for (auto _ : state)
        for (auto && query : reads)
            benchmark::DoNotOptimize(agent.bulk_count(query | hash_adaptor));
}

BENCHMARK(bulk_count)->RangeMultiplier(2)->Range(64, 65536);

BENCHMARK_MAIN();

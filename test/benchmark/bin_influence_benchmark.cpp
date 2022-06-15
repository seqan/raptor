// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <filesystem>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/algorithm/detail/execution_handler_parallel.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/utility/views/slice.hpp>
#include <seqan3/utility/views/zip.hpp>

static constexpr size_t operator""_MiB(unsigned long long int number)
{
    return number << 23;
}

static constexpr size_t operator""_GiB(unsigned long long int number)
{
    return number << 33;
}

#if 1
static constexpr size_t const genome_size{5000};
static constexpr size_t const read_size{250};
static constexpr size_t const read_count{1000};
static constexpr size_t const max_ibf_size{1_MiB};
#else
static constexpr size_t const genome_size{4'300'000'000};
static constexpr size_t const read_size{250};
static constexpr size_t const read_count{1ULL << 20};
static constexpr size_t const max_ibf_size{10_GiB};
#endif

static constexpr size_t const hash_num{4u};
static constexpr size_t const window_size{24};
static constexpr size_t const kmer_size{20};
static std::filesystem::path base_path{"/dev/shm/seiler/bin_influence"};
static constexpr size_t const construct_threads{32}; // Only applies to construction

static std::vector<seqan3::dna4> const genome{seqan3::test::generate_sequence<seqan3::dna4>(genome_size, 0, 0)};
static std::vector<std::vector<seqan3::dna4>> const reads{
    [](auto const & genome)
    {
        std::vector<std::vector<seqan3::dna4>> result(read_count);
        size_t i{};
        for (auto && read_start :
             seqan3::test::generate_numeric_sequence<size_t>(read_count, 0, genome_size - read_size + 1, 0))
        {
            auto v = genome | seqan3::views::slice(read_start, read_start + read_size);
            result[i++].assign(v.begin(), v.end());
        }
        return result;
    }(genome)};

using ibf_t = seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>;

static constexpr uint64_t adjust_seed(uint8_t const kmer_size, uint64_t const seed = 0x8F3F73B5CF1C9ADEULL) noexcept
{
    return seed >> (64u - 2u * kmer_size);
}

static constexpr size_t compute_bin_size(size_t const max_bin_size, double const fpr)
{
    double const numerator{-static_cast<double>(max_bin_size * hash_num)};
    double const denominator{std::log(1 - std::exp(std::log(fpr) / hash_num))};
    double const result{std::ceil(numerator / denominator)};
    return result;
}

static std::vector<size_t> cardinality(size_t const bin_count, auto && hash_adaptor)
{
    std::vector<size_t> cardinalities(bin_count);

    size_t const chunked_genome_size{(genome_size + bin_count - 1) / bin_count};
    auto chunked_genomes = genome | seqan3::views::chunk(chunked_genome_size);

    size_t const workload_size = std::clamp<size_t>(std::bit_ceil(bin_count / construct_threads), 8u, 64u);
    auto workload = seqan3::views::zip(chunked_genomes, std::views::iota(0u)) | seqan3::views::chunk(workload_size);

    auto worker = [&cardinalities, &hash_adaptor](auto && payload, auto &&)
    {
        for (auto && [sequence, bin_number] : payload)
            cardinalities[bin_number] += std::ranges::distance(sequence | hash_adaptor);
    };

    seqan3::detail::execution_handler_parallel executioner{construct_threads};
    executioner.bulk_execute(std::move(worker), std::move(workload), []() {});

    return cardinalities;
}

static ibf_t construct_ibf(size_t const bin_count, auto && hash_adaptor, double const fpr)
{
    size_t const bin_size{compute_bin_size(std::ranges::max(cardinality(bin_count, hash_adaptor)), fpr)};

    if (bin_size * bin_count > max_ibf_size)
        throw std::runtime_error{"Resulting IBF would be too big. " + std::to_string(bin_size * bin_count)};

    ibf_t ibf{seqan3::bin_count{bin_count}, seqan3::bin_size{bin_size}, seqan3::hash_function_count{hash_num}};

    size_t const chunked_genome_size{(genome_size + bin_count - 1) / bin_count};
    auto chunked_genomes = genome | seqan3::views::chunk(chunked_genome_size);

    size_t const workload_size = std::clamp<size_t>(std::bit_ceil(bin_count / construct_threads), 8u, 64u);
    auto workload = seqan3::views::zip(chunked_genomes, std::views::iota(0u)) | seqan3::views::chunk(workload_size);

    auto worker = [&ibf, &hash_adaptor](auto && payload, auto &&)
    {
        for (auto && [sequence, bin_number] : payload)
            for (auto && hash : sequence | hash_adaptor)
                ibf.emplace(hash, seqan3::bin_index{bin_number});
    };

    seqan3::detail::execution_handler_parallel executioner{construct_threads};
    executioner.bulk_execute(std::move(worker), std::move(workload), []() {});

    return ibf;
}

static void bulk_count(benchmark::State & state, double && fpr)
{
    size_t const bin_count = static_cast<size_t>(state.range(0));
    auto hash_adaptor = seqan3::views::minimiser_hash(seqan3::ungapped{kmer_size},
                                                      seqan3::window_size{window_size},
                                                      seqan3::seed{adjust_seed(kmer_size)});

    std::filesystem::path ibf_path = base_path;
    ibf_path /= std::to_string(window_size) + "_" + std::to_string(kmer_size) + "_" + std::to_string(fpr) + "_"
              + std::to_string(bin_count) + ".ibf";

    ibf_t ibf{};
    if (std::filesystem::exists(ibf_path))
    {
        std::ifstream is{ibf_path, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        iarchive(ibf);
    }
    else
    {
        ibf = std::move(construct_ibf(bin_count, hash_adaptor, fpr));
        std::ofstream os{ibf_path, std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(ibf);
    }

    auto agent = ibf.counting_agent<uint16_t>();
    for (auto _ : state)
        for (auto && query : reads)
            benchmark::DoNotOptimize(agent.bulk_count(query | hash_adaptor));
}

BENCHMARK_CAPTURE(bulk_count, "0.0001", 0.0001)->RangeMultiplier(2)->Range(64, 32768);
// BENCHMARK_CAPTURE(bulk_count, "0.0005", 0.0005)->RangeMultiplier(2)->Range(64, 32768);
// BENCHMARK_CAPTURE(bulk_count, "0.0025", 0.0025)->RangeMultiplier(2)->Range(64, 32768);
BENCHMARK_CAPTURE(bulk_count, "0.0125", 0.0125)->RangeMultiplier(2)->Range(64, 32768);
// BENCHMARK_CAPTURE(bulk_count, "0.0625", 0.0625)->RangeMultiplier(2)->Range(64, 32768);
BENCHMARK_CAPTURE(bulk_count, "0.3125", 0.3125)->RangeMultiplier(2)->Range(64, 32768);

BENCHMARK_MAIN();

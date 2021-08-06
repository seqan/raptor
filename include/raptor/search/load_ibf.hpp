#pragma once

#include <chrono>
#include <seqan3/std/filesystem>

#include <raptor/shared.hpp>

namespace raptor
{

template <typename t>
void load_ibf(t & ibf, search_arguments const & arguments, size_t const part, double & ibf_io_time)
{
    static uint8_t kmer_size{};
    static uint32_t window_size{};
    static std::vector<std::vector<std::string>> bin_path{};

    std::filesystem::path ibf_file{arguments.ibf_file};
    ibf_file += "_" + std::to_string(part);

    std::ifstream is{ibf_file, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};

    auto start = std::chrono::high_resolution_clock::now();
    iarchive(kmer_size);
    iarchive(window_size);
    iarchive(bin_path);
    iarchive(ibf);
    auto end = std::chrono::high_resolution_clock::now();

    ibf_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
}

template <typename t>
void load_ibf(t & ibf, search_arguments const & arguments, double & ibf_io_time)
{
    static uint8_t kmer_size{};
    static uint32_t window_size{};
    static std::vector<std::vector<std::string>> bin_path{};

    std::ifstream is{arguments.ibf_file, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};

    auto start = std::chrono::high_resolution_clock::now();
    iarchive(kmer_size);
    iarchive(window_size);
    iarchive(bin_path);
    iarchive(ibf);
    auto end = std::chrono::high_resolution_clock::now();

    ibf_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
}

} // namespace raptor

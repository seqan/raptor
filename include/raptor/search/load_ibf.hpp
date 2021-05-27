#pragma once

#include <chrono>
#include <seqan3/std/filesystem>

#include <raptor/shared.hpp>

namespace raptor
{

template <typename t>
void load_ibf(t & ibf, search_arguments const & arguments, size_t const part, double & ibf_io_time)
{
    std::filesystem::path ibf_file{arguments.ibf_file};
    ibf_file += "_" + std::to_string(part);
    std::ifstream is{ibf_file, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    auto start = std::chrono::high_resolution_clock::now();
    iarchive(ibf);
    auto end = std::chrono::high_resolution_clock::now();
    ibf_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
}

} // namespace raptor

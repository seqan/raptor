#pragma once

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

#include <raptor/shared.hpp>

namespace raptor
{

template <auto layout>
static inline void store_index(std::filesystem::path const & path,
                               seqan3::interleaved_bloom_filter<layout> const & ibf,
                               build_arguments const & arguments)
{
    std::ofstream os{path, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(arguments.kmer_size);
    oarchive(arguments.window_size);
    oarchive(ibf);
}

} // namespace raptor


// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::build_ibf.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <hibf/build/bin_size_in_bits.hpp>

#include <raptor/build/index_factory.hpp>
#include <raptor/build/max_count_per_partition.hpp>
#include <raptor/build/partition_config.hpp>
#include <raptor/build/store_index.hpp>

namespace raptor
{

void build_ibf(build_arguments const & arguments)
{
    if (arguments.parts == 1u)
    {
        index_factory factory{arguments};
        auto index = factory();
        arguments.store_index_timer.start();
        store_index(arguments.out_path, std::move(index));
        arguments.store_index_timer.stop();
    }
    else
    {
        partition_config const cfg{arguments.parts};
        index_factory factory{arguments, cfg};
        std::vector<size_t> const kmers_per_partition = max_count_per_partition(cfg, arguments);

        for (size_t part = 0; part < arguments.parts; ++part)
        {
            arguments.bits = seqan::hibf::build::bin_size_in_bits(
                {.fpr = arguments.fpr, .hash_count = arguments.hash, .elements = kmers_per_partition[part]});
            auto index = factory(part);
            std::filesystem::path out_path{arguments.out_path};
#if HIBF_WORKAROUND_GCC_BOGUS_MEMCPY
#    pragma GCC diagnostic push
#    pragma GCC diagnostic ignored "-Wrestrict"
#endif // HIBF_WORKAROUND_GCC_BOGUS_MEMCPY
            out_path += "_" + std::to_string(part);
#if HIBF_WORKAROUND_GCC_BOGUS_MEMCPY
#    pragma GCC diagnostic pop
#endif // HIBF_WORKAROUND_GCC_BOGUS_MEMCPY
            arguments.store_index_timer.start();
            store_index(out_path, std::move(index));
            arguments.store_index_timer.stop();
        }
    }
}

} // namespace raptor

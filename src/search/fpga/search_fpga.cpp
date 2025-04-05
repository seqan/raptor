// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-FileCopyrightText: 2020-2025 Thomas Steinke & Zuse Institute Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <cmath>
#include <fstream>
#include <vector>

#include <hibf/misc/unreachable.hpp>

#include <raptor/search/fpga/min_ibf_fpga_oneapi.hpp>
#include <raptor/search/fpga/search_fpga.hpp>
#include <raptor/search/sync_out.hpp>
#include <raptor/threshold/threshold.hpp>

#ifndef FPGA_BINS
#    error "FPGA_BINS is not defined."
#endif

namespace raptor
{

class fpga_thresholder : raptor::threshold::threshold
{
private:
    using base_t = raptor::threshold::threshold;

public:
    using base_t::base_t;

    std::pair<size_t, size_t> get_minmax() const
    {
        return std::make_pair(minimal_number_of_minimizers, maximal_number_of_minimizers);
    }

    std::vector<size_t> get_thresholds() const
    {
        assert(threshold_kind == threshold_kinds::probabilistic);

        std::vector<size_t> result;

        std::ranges::transform(precomp_thresholds,
                               precomp_correction,
                               std::back_inserter(result),
                               [](size_t const threshold, size_t const correction)
                               {
                                   return std::max<size_t>(1u, threshold + correction);
                               });

        return result;
    }
};

void search_fpga(search_arguments const & arguments)
{
    std::vector<size_t> thresholds;
    size_t minimal_number_of_minimizers{};
    size_t maximal_number_of_minimizers{};
    {
        fpga_thresholder thresholder{arguments.make_threshold_parameters()};
        thresholds = thresholder.get_thresholds();
        std::tie(minimal_number_of_minimizers, maximal_number_of_minimizers) = thresholder.get_minmax();
    }

    size_t const bins{arguments.bin_path.size()};
    size_t const technical_bins{seqan::hibf::next_multiple_of_64(bins)};
    assert(bins == technical_bins); // Todo: Important?

    constexpr bool profile = true;

    // size_t const chunk_bits = std::min<size_t>(technical_bins, MAX_BUS_WIDTH);
    // assert(MAX_BUS_WIDTH <= 512);

    auto process = [&]<size_t technical_bins>()
    {
        min_ibf_fpga_oneapi<technical_bins, profile> ibf(arguments.index_file,
                                                         minimal_number_of_minimizers,
                                                         maximal_number_of_minimizers,
                                                         std::move(thresholds),
                                                         arguments.buffer,
                                                         arguments.kernels);

        {
            // Todo: Not so elegant. min_ibf_fpga_oneapi also needs to open the file in append mode for this to work.
            sync_out synced_out{arguments};
            synced_out.write_header(arguments, ibf.hash_function_count());
        }

        ibf.count(arguments.query_file, arguments.out_file);
    };

    constexpr auto allowed_bins = std::to_array<size_t>({FPGA_BINS});
    static_assert(!std::empty(allowed_bins));

    // Same as `for (size_t bin : allowed_bins) { if (technical_bins == bin) { process.template operator()<bin>(); break; } }`
    // but constexpr evaluated, such that templates are instantiated.
    [&]<size_t... idx>(std::index_sequence<idx...>)
    {
        ((technical_bins == allowed_bins[idx] ? (process.template operator()<allowed_bins[idx]>(), void()) : void()),
         ...);
    }(std::make_index_sequence<allowed_bins.size()>{});
}

} // namespace raptor

// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides raptor::raptor_upgrade.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <raptor/argument_parsing/compute_bin_size.hpp>
#include <raptor/argument_parsing/validators.hpp>
#include <raptor/build/max_count_per_partition.hpp>
#include <raptor/upgrade/index_upgrader.hpp>
#include <raptor/upgrade/upgrade.hpp>

namespace raptor
{

void raptor_upgrade(upgrade_arguments & arguments)
{
    std::variant<index_upgrader<seqan3::data_layout::uncompressed>, index_upgrader<seqan3::data_layout::compressed>>
        upgrader{};

    if (arguments.parts == 1u)
    {
        size_t const max_count = std::isnan(arguments.fpr) ? max_bin_count(arguments) : 0u;

        if (arguments.compressed)
            upgrader = index_upgrader<seqan3::data_layout::compressed>{arguments, max_count};
        else
            upgrader = index_upgrader<seqan3::data_layout::uncompressed>{arguments, max_count};

        std::visit(
            [](auto & u)
            {
                u.upgrade();
            },
            upgrader);
    }
    else
    {
        partition_config const cfg{arguments.parts};
        std::vector<size_t> count_per_partition =
            std::isnan(arguments.fpr) ? max_count_per_partition(cfg, arguments) : std::vector<size_t>{};
        std::string const index_path_base = arguments.index_file.string() + '_';
        std::string const output_path_base = arguments.output_file.string() + '_';

        for (size_t part{0}; part < arguments.parts; ++part)
        {
            arguments.index_file = index_path_base + std::to_string(part);
            arguments.output_file = output_path_base + std::to_string(part);

            size_t const max_count = std::isnan(arguments.fpr) ? count_per_partition[part] : 0u;

            if (arguments.compressed)
                upgrader = index_upgrader<seqan3::data_layout::compressed>{arguments, max_count}; // GCOVR_EXCL_LINE
            else
                upgrader = index_upgrader<seqan3::data_layout::uncompressed>{arguments, max_count};

            std::visit(
                [](auto & u)
                {
                    u.upgrade();
                },
                upgrader);
        }
    }
}

} // namespace raptor

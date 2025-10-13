// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::raptor_upgrade.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <raptor/argument_parsing/upgrade_arguments.hpp> // for upgrade_arguments
#include <raptor/upgrade/upgrade.hpp>                    // for raptor_upgrade

namespace raptor
{

void raptor_upgrade(upgrade_arguments & /* arguments */) // GCOVR_EXCL_LINE
{
    // if (arguments.parts == 1u)
    // {
    //     size_t const max_count = std::isnan(arguments.fpr) ? max_bin_count(arguments) : 0u;

    //     index_upgrader upgrader{arguments, max_count};
    //     upgrader.upgrade();
    // }
    // else
    // {
    //     partition_config const cfg{arguments.parts};
    //     std::vector<size_t> count_per_partition =
    //         std::isnan(arguments.fpr) ? max_count_per_partition(cfg, arguments) : std::vector<size_t>{};
    //     std::string const index_path_base = arguments.index_file.string() + '_';
    //     std::string const output_path_base = arguments.output_file.string() + '_';

    //     for (size_t part{0}; part < arguments.parts; ++part)
    //     {
    //         arguments.index_file = index_path_base + std::to_string(part);
    //         arguments.output_file = output_path_base + std::to_string(part);

    //         size_t const max_count = std::isnan(arguments.fpr) ? count_per_partition[part] : 0u;

    //         index_upgrader upgrader{arguments, max_count};
    //         upgrader.upgrade();
    //     }
    // }
}

} // namespace raptor

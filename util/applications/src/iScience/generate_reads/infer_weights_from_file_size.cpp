// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "cmd_arguments.hpp"

namespace raptor::util::generate_reads
{

void infer_weights_from_file_size(cmd_arguments & arguments)
{
#pragma omp parallel for schedule(static) num_threads(arguments.threads)
    for (size_t i = 0; i < arguments.number_of_bins; ++i)
    {
        arguments.number_of_reads_per_bin[i] = std::filesystem::file_size(arguments.bin_path[i]);
    }
}

} // namespace raptor::util::generate_reads

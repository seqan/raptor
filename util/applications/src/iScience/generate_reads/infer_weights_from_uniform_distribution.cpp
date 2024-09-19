// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <algorithm>

#include "cmd_arguments.hpp"

namespace raptor::util::generate_reads
{

void infer_weights_from_uniform_distribution(cmd_arguments & arguments)
{
    std::ranges::fill(arguments.number_of_reads_per_bin, arguments.number_of_reads);
}

} // namespace raptor::util::generate_reads

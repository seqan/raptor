// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <hibf/sketch/compute_sketches.hpp>
#include <hibf/sketch/hyperloglog.hpp>

#include <raptor/file_reader.hpp>

#include "cmd_arguments.hpp"

namespace raptor::util::generate_reads
{

void infer_weights_from_hll_sketches(cmd_arguments & arguments)
{
    auto input_fn = [&arguments](size_t const ub_id, seqan::hibf::insert_iterator it)
    {
        static raptor::file_reader<raptor::file_types::sequence> const reader{seqan3::ungapped{32u}, 32u};
        reader.hash_into(arguments.bin_path[ub_id], it);
    };

    seqan::hibf::config config{.input_fn = input_fn,
                               .number_of_user_bins = arguments.number_of_bins,
                               .threads = arguments.threads,
                               .sketch_bits = 12};

    std::vector<seqan::hibf::sketch::hyperloglog> hyperloglog_sketches;
    seqan::hibf::sketch::compute_sketches(config, hyperloglog_sketches);

    assert(hyperloglog_sketches.size() == arguments.number_of_bins);

#pragma omp parallel for schedule(static) num_threads(arguments.threads)
    for (size_t i = 0; i < arguments.number_of_bins; ++i)
    {
        arguments.number_of_reads_per_bin[i] = hyperloglog_sketches[i].estimate();
    }
}

} // namespace raptor::util::generate_reads

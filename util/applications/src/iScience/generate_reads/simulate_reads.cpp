// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <random>

#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>

#include <hibf/misc/divide_and_ceil.hpp>

#include "cmd_arguments.hpp"
#include "seqan3_sequence_io.hpp"

namespace raptor::util::generate_reads
{

size_t simulate_reads(cmd_arguments const & arguments)
{
    using input_file_t = seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>>;
    std::vector<seqan3::phred42> const quality(arguments.read_length, seqan3::assign_rank_to(40u, seqan3::phred42{}));
    seqan3::sequence_file_output fout{arguments.output_filename};

    std::uniform_int_distribution<size_t> read_error_position_dis(0u, arguments.read_length - 1u);
    std::uniform_int_distribution<uint8_t> dna4_rank_dis(0u, 3u);

    size_t total_reads{};

#pragma omp parallel for schedule(dynamic) num_threads(arguments.threads)
    for (size_t bin_number = 0; bin_number < arguments.bin_path.size(); ++bin_number)
    {
        std::mt19937_64 rng{bin_number};

        auto const & bin_file = arguments.bin_path[bin_number];

        size_t const number_of_records{count_records_in_fasta(bin_file)};
        size_t const reads_per_record =
            seqan::hibf::divide_and_ceil(arguments.number_of_reads_per_bin[bin_number], number_of_records);

        std::vector<std::vector<seqan3::dna4>> reads(reads_per_record * number_of_records);

        size_t read_counter{};

        for (auto const & [seq] : input_file_t{bin_file})
        {
            size_t const reference_length = std::ranges::size(seq);
            if (reference_length < arguments.read_length)
                continue;
            std::uniform_int_distribution<size_t> read_start_dis(0u, reference_length - arguments.read_length);

            for (size_t read_number = 0u; read_number < reads_per_record; ++read_number, ++read_counter)
            {
                auto & read = reads[read_counter];
                size_t const read_start_pos = read_start_dis(rng);
                auto read_span = std::span{seq.data() + read_start_pos, arguments.read_length};
                read.assign(read_span.begin(), read_span.end());

                for (uint8_t error_count = 0; error_count < arguments.errors; ++error_count)
                {
                    size_t const error_pos = read_error_position_dis(rng);
                    seqan3::dna4 const current_base = read[error_pos];
                    seqan3::dna4 new_base = current_base;
                    while (new_base == current_base)
                        seqan3::assign_rank_to(dna4_rank_dis(rng), new_base);
                    read[error_pos] = new_base;
                }
            }
        }

        std::string const bin_filename = bin_file.filename();

#pragma omp critical
        {
            total_reads += read_counter;
            for (size_t read_number = 0u; read_number < read_counter; ++read_number)
            {
                fout.emplace_back(reads[read_number], bin_filename + "_" + std::to_string(read_number), quality);
            }
        }
    }

    return total_reads;
}

} // namespace raptor::util::generate_reads

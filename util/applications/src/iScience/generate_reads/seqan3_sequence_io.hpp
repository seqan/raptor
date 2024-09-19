// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>

namespace raptor::util::generate_reads
{

struct char_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = char;
    using sequence_legal_alphabet = char;
};

struct dna4_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

size_t count_records_in_fasta(std::filesystem::path const & filename)
{
    seqan3::sequence_file_input<char_traits, seqan3::fields<seqan3::field::id>> fin{filename};
    return std::ranges::distance(fin);
}

} // namespace raptor::util::generate_reads

// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <raptor/argument_parsing/top_level.hpp>

namespace raptor
{

void init_top_level_parser(seqan3::argument_parser & parser)
{
    init_shared_meta(parser);
    parser.info.description.emplace_back("Raptor is a system for approximately searching many queries such as "
                                         "next-generation sequencing reads or transcripts in large collections of "
                                         "nucleotide sequences. Raptor uses winnowing minimizers to define a set of "
                                         "representative k-mers, an extension of the interleaved Bloom filters (IBFs) "
                                         "as a set membership data structure and probabilistic thresholding for "
                                         "minimizers. Our approach allows compression and partitioning of the IBF to "
                                         "enable the effective use of secondary memory.");
};

} // namespace raptor

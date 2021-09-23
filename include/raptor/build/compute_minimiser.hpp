// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <raptor/argument_parsing/build_arguments.hpp>

namespace raptor
{

bool check_for_fasta_format(std::vector<std::string> const & valid_extensions, std::string const & file_path);

void compute_minimiser(build_arguments const & arguments);

} // namespace raptor

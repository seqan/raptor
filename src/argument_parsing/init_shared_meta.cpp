// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <raptor/argument_parsing/init_shared_meta.hpp>

namespace raptor
{

void init_shared_meta(sharg::parser & parser)
{
    parser.info.app_name = "Raptor";
    parser.info.author = "Enrico Seiler";
    parser.info.citation = "Raptor: A fast and space-efficient pre-filter for querying very large collections of "
                           "nucleotide sequences; Enrico Seiler, Svenja Mehringer, Mitra Darvish, Etienne Turc, "
                           "and Knut Reinert; iScience 2021 24 (7): 102782. doi: "
                           "https://doi.org/10.1016/j.isci.2021.102782";
    parser.info.date = RAPTOR_DATE;
    parser.info.email = "enrico.seiler@fu-berlin.de";
    parser.info.long_copyright = R"(BSD 3-Clause License

Copyright (c) 2022, Enrico Seiler
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.)";
    parser.info.short_copyright = "BSD 3-Clause License";
    parser.info.short_description =
        "A fast and space-efficient pre-filter for querying very large collections of nucleotide sequences.";
    parser.info.url = "https://github.com/seqan/raptor";
    parser.info.version = RAPTOR_VERSION;
}

} // namespace raptor

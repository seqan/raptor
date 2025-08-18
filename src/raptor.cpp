// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <raptor/argument_parsing/build_parsing.hpp>
#include <raptor/argument_parsing/prepare_parsing.hpp>
#include <raptor/argument_parsing/search_parsing.hpp>
#include <raptor/argument_parsing/update_parsing.hpp>
#include <raptor/argument_parsing/upgrade_parsing.hpp>
#include <raptor/layout/raptor_layout.hpp>
#include <raptor/raptor.hpp>

inline void set_metadata(sharg::parser_meta_data & info)
{
    info.version = RAPTOR_VERSION;
    info.short_description =
        "A fast and space-efficient pre-filter for querying very large collections of nucleotide sequences";
    info.author = "Enrico Seiler";
    info.email = "enrico.seiler@fu-berlin.de";
    info.date = RAPTOR_DATE;
    info.url = "https://github.com/seqan/raptor";
    info.short_copyright = "BSD 3-Clause License";
    info.long_copyright = R"(BSD 3-Clause License

Copyright (c) 2006-2025, Knut Reinert & Freie Universität Berlin
Copyright (c) 2016-2025, Knut Reinert & MPI für molekulare Genetik
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
    info.citation = {"Raptor: A fast and space-efficient pre-filter for querying very large collections of "
                     "nucleotide sequences; Enrico Seiler, Svenja Mehringer, Mitra Darvish, Etienne Turc, "
                     "and Knut Reinert; iScience 2021 24 (7): 102782. doi: "
                     "https://doi.org/10.1016/j.isci.2021.102782",
                     "Hierarchical Interleaved Bloom Filter: "
                     "enabling ultrafast, approximate sequence queries; Svenja Mehringer, Enrico Seiler, Felix "
                     "Droop, Mitra Darvish, René Rahn, Martin Vingron, and Knut Reinert; Genome Biol 24, 131 "
                     "(2023). doi: https://doi.org/10.1186/s13059-023-02971-4"};
    info.description.emplace_back("Raptor is a system for approximately searching many queries such as "
                                  "next-generation sequencing reads or transcripts in large collections of "
                                  "nucleotide sequences. Raptor uses winnowing minimizers to define a set of "
                                  "representative k-mers, an extension of the interleaved Bloom filters (IBFs) "
                                  "as a set membership data structure and probabilistic thresholding for "
                                  "minimizers. Our approach allows compression and partitioning of the IBF to "
                                  "enable the effective use of secondary memory.");
}

int main(int argc, char ** argv)
{
    try
    {
        sharg::parser top_level_parser{"Raptor",
                                       argc,
                                       argv,
                                       sharg::update_notifications::on,
                                       {"build", "layout", "prepare", "search", "update", "upgrade"}};
        set_metadata(top_level_parser.info);

        top_level_parser.parse();

        sharg::parser & sub_parser = top_level_parser.get_sub_parser();
        if (sub_parser.info.app_name == std::string_view{"Raptor-build"})
            raptor::build_parsing(sub_parser);
        if (sub_parser.info.app_name == std::string_view{"Raptor-layout"})
            raptor::chopper_layout(sub_parser);
        if (sub_parser.info.app_name == std::string_view{"Raptor-prepare"})
            raptor::prepare_parsing(sub_parser);
        if (sub_parser.info.app_name == std::string_view{"Raptor-search"})
            raptor::search_parsing(sub_parser);
        if (sub_parser.info.app_name == std::string_view{"Raptor-update"})
            raptor::update_parsing(sub_parser);
        if (sub_parser.info.app_name == std::string_view{"Raptor-upgrade"})
            raptor::upgrade_parsing(sub_parser); // GCOVR_EXCL_LINE
    }
    catch (std::exception const & ext)
    {
        std::cerr << "[Error] " << ext.what() << '\n';
        std::exit(-1);
    }

    return 0;
}

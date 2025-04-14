// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <cstdlib>     // for exit
#include <exception>   // for exception
#include <iostream>    // for basic_ostream, operator<<, basic_ios, cerr
#include <string>      // for basic_string, char_traits
#include <string_view> // for basic_string_view, operator==, string_view
#include <vector>      // for vector

#include <sharg/auxiliary.hpp> // for parser_meta_data, update_notifications
#include <sharg/parser.hpp>    // for parser

#include <raptor/argument_parsing/build_parsing.hpp>    // for build_parsing
#include <raptor/argument_parsing/init_shared_meta.hpp> // for init_shared_meta
#include <raptor/argument_parsing/prepare_parsing.hpp>  // for prepare_parsing
#include <raptor/argument_parsing/search_parsing.hpp>   // for search_parsing
#include <raptor/argument_parsing/update_parsing.hpp>   // for update_parsing
#include <raptor/argument_parsing/upgrade_parsing.hpp>  // for upgrade_parsing
#include <raptor/layout/raptor_layout.hpp>              // for chopper_layout

int main(int argc, char ** argv)
{
    try
    {
        sharg::parser top_level_parser{"Raptor",
                                       argc,
                                       argv,
                                       sharg::update_notifications::on,
                                       {"build", "layout", "prepare", "search", "update", "upgrade"}};
        raptor::init_shared_meta(top_level_parser);
        top_level_parser.info.description.emplace_back(
            "Raptor is a system for approximately searching many queries such as "
            "next-generation sequencing reads or transcripts in large collections of "
            "nucleotide sequences. Raptor uses winnowing minimizers to define a set of "
            "representative k-mers, an extension of the interleaved Bloom filters (IBFs) "
            "as a set membership data structure and probabilistic thresholding for "
            "minimizers. Our approach allows compression and partitioning of the IBF to "
            "enable the effective use of secondary memory.");

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

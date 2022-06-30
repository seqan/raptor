// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <raptor/argument_parsing/build_parsing.hpp>
#include <raptor/argument_parsing/init_shared_meta.hpp>
#include <raptor/argument_parsing/search_parsing.hpp>
#include <raptor/argument_parsing/upgrade_parsing.hpp>
#include <raptor/raptor.hpp>

int main(int argc, char ** argv)
{
    try
    {
        sharg::parser top_level_parser{"raptor",
                                       argc,
                                       argv,
                                       sharg::update_notifications::on,
                                       {"build", "search", "socks", "upgrade"}};
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
        if (sub_parser.info.app_name == std::string_view{"raptor-build"})
            raptor::build_parsing(sub_parser, false);
        if (sub_parser.info.app_name == std::string_view{"raptor-search"})
            raptor::search_parsing(sub_parser, false);
        if (sub_parser.info.app_name == std::string_view{"raptor-socks"})
        {
            sharg::parser socks_parser{"socks",
                                       argc - 1,
                                       argv + 1,
                                       sharg::update_notifications::off,
                                       {"build", "lookup-kmer"}};
            socks_parser.parse();
            sharg::parser & socks_sub_parser = socks_parser.get_sub_parser();
            if (socks_sub_parser.info.app_name == std::string_view{"socks-build"})
                raptor::build_parsing(socks_sub_parser, true);
            if (socks_sub_parser.info.app_name == std::string_view{"socks-lookup-kmer"})
                raptor::search_parsing(socks_sub_parser, true);
        }
        if (sub_parser.info.app_name == std::string_view{"raptor-upgrade"})
            raptor::upgrade_parsing(sub_parser);
    }
    catch (sharg::parser_error const & ext)
    {
        std::cerr << "[Error] " << ext.what() << '\n';
        std::exit(-1);
    }

    return 0;
}

// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <seqan3/argument_parser/all.hpp>

#include <raptor/argument_parsing/build.hpp>
#include <raptor/argument_parsing/search.hpp>
#include <raptor/argument_parsing/shared.hpp>
#include <raptor/argument_parsing/upgrade.hpp>
#include <raptor/argument_parsing/top_level.hpp>
#include <raptor/raptor.hpp>

int main(int argc, char ** argv)
{
    try
    {
        seqan3::argument_parser top_level_parser{"raptor", argc, argv, seqan3::update_notifications::on, {"build", "search", "socks", "upgrade"}};
        raptor::init_top_level_parser(top_level_parser);

        raptor::try_parsing(top_level_parser);

        seqan3::argument_parser & sub_parser = top_level_parser.get_sub_parser();
        if (sub_parser.info.app_name == std::string_view{"raptor-build"})
            raptor::run_build(sub_parser, false);
        if (sub_parser.info.app_name == std::string_view{"raptor-search"})
            raptor::run_search(sub_parser, false);
        if (sub_parser.info.app_name == std::string_view{"raptor-socks"})
        {
            seqan3::argument_parser socks_parser{"socks", argc - 1, argv + 1, seqan3::update_notifications::off, {"build", "lookup-kmer"}};
            raptor::try_parsing(socks_parser);
            seqan3::argument_parser & socks_sub_parser = socks_parser.get_sub_parser();
            if (socks_sub_parser.info.app_name == std::string_view{"socks-build"})
                raptor::run_build(socks_sub_parser, true);
            if (socks_sub_parser.info.app_name == std::string_view{"socks-lookup-kmer"})
                raptor::run_search(socks_sub_parser, true);
        }
        if (sub_parser.info.app_name == std::string_view{"raptor-upgrade"})
            raptor::run_upgrade(sub_parser);
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        std::cerr << "[Error] " << ext.what() << '\n';
        std::exit(-1);
    }

    return 0;
}

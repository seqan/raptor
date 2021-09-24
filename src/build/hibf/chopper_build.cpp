// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <raptor/build/hibf/build_config.hpp>
#include <raptor/build/hibf/build_data.hpp>
#include <raptor/build/hibf/chopper_build.hpp>
#include <raptor/build/hibf/create_ibfs_from_chopper_pack.hpp>

namespace raptor::hibf
{

void initialize_argument_parser(seqan3::argument_parser & parser, build_config & config)
{
    parser.add_option(config.chopper_pack_filename, 'p', "pack-file", "Provide the file produced by chopper pack.");
    parser.add_option(config.k, 'k', "kmer-size", "The kmer size to build kmers.");
    parser.add_option(config.hash_funs, '\0', "hash-functions", "The number of hash functions to use for the IBF.");
    parser.add_option(config.overlap, 'l', "overlap", "The overlap between split regions of the same sequence.");
    parser.add_option(config.FPR, 'r', "false-positive-rate", "The minimum false positive rate of every IBF.");
    parser.add_option(config.threads, 't', "threads", "The number of threads to use.");
    parser.add_option(config.output_filename, 'o', "out-prefix", "Output filename.");
    parser.add_flag(config.verbose, 'v', "verbose", "Output logging/progress information.");
}

void chopper_build(seqan3::argument_parser & parser)
{
    build_config config{};
    initialize_argument_parser(parser, config);

    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        std::cerr << "[CHOPPER BUILD ERROR] " << ext.what() << "\n";
    }

    build_data data{};
    create_ibfs_from_chopper_pack(data, config);

    // write vector of ibfs to file
    {
        std::ofstream fout(config.output_filename, std::ios::binary);

        if (!fout.good() || !fout.is_open())
            throw std::runtime_error{"Could not open " + config.output_filename + " for writing."};

        cereal::BinaryOutputArchive archive(fout);
        archive(data.hibf);
        archive(config.k);
    }
}

} // namespace raptor::hibf



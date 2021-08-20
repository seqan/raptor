// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <raptor/build/compute_minimiser.hpp>

namespace raptor
{

bool check_for_fasta_format(std::vector<std::string> const & valid_extensions, std::string const & file_path)
{

    auto case_insensitive_string_ends_with = [&] (std::string_view str, std::string_view suffix)
    {
        size_t const suffix_length{suffix.size()};
        size_t const str_length{str.size()};
        return suffix_length > str_length ?
               false :
               std::ranges::equal(str.substr(str_length - suffix_length), suffix, [] (char const chr1, char const chr2)
               {
                   return std::tolower(chr1) == std::tolower(chr2);
               });
    };

    auto case_insensitive_ends_with = [&] (std::string const & ext)
    {
        return case_insensitive_string_ends_with(file_path, ext);
    };

    return std::ranges::find_if(valid_extensions, case_insensitive_ends_with) != valid_extensions.end();
}

void compute_minimiser(build_arguments const & arguments)
{
    auto minimiser_view = seqan3::views::minimiser_hash(seqan3::ungapped{arguments.kmer_size},
                                                        seqan3::window_size{arguments.window_size},
                                                        seqan3::seed{adjust_seed(arguments.kmer_size)});

    uint16_t const default_cutoff{50};

    // Cutoffs and bounds from Mantis
    // Mantis ignores k-mers which appear less than a certain cutoff. The cutoff is based on the file size of a
    // gzipped fastq file. Small files have only a cutoff of 1 while big files have a cutoff value of 50.
    std::array<uint16_t, 4> const cutoffs{1, 3, 10, 20};
    std::array<uint64_t, 4> const cutoff_bounds{314'572'800, 524'288'000, 1'073'741'824, 3'221'225'472};

    auto worker = [&] (auto && zipped_view, auto &&)
    {
        robin_hood::unordered_map<uint64_t, uint8_t> minimiser_table{};
        uint64_t count{0};
        uint16_t cutoff{0};

        for (auto && [file_names, bin_number] : zipped_view)
        {
            for (auto && file_name : file_names)
            {
                seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>> fin{file_name};

                for (auto & [seq] : fin)
                    for (auto && hash : seq | minimiser_view)
                        minimiser_table[hash] = std::min<uint8_t>(254u, minimiser_table[hash] + 1);
                        // The hash table stores how often a minimiser appears. It does not matter whether a minimiser appears
                        // 50 times or 2000 times, it is stored regardless because the biggest cutoff value is 50. Hence,
                        // the hash table stores only values up to 254 to save memory.
            }

            std::filesystem::path const file_name{file_names[0]};

            if (!arguments.disable_cutoffs)
            {
                // Since the curoffs are based on the filesize of a gzipped fastq file, we try account for the other cases:
                // We multiply by two if we have fasta input.
                // We divide by 3 if the input is not compressed.
                bool const is_compressed = file_name.extension() == ".gz" || file_name.extension() == ".bgzf" || file_name.extension() == ".bz2";
                bool const is_fasta = is_compressed ? check_for_fasta_format(seqan3::format_fasta::file_extensions, file_name.stem())
                                                    : check_for_fasta_format(seqan3::format_fasta::file_extensions, file_name.extension());
                size_t const filesize = std::filesystem::file_size(file_name) * (is_fasta ? 2 : 1) / (is_compressed ? 1 : 3);

                cutoff = default_cutoff;
                for (size_t k = 0; k < cutoff_bounds.size(); ++k)
                {
                    if (filesize <= cutoff_bounds[k])
                    {
                        cutoff = cutoffs[k];
                        break;
                    }
                }
            }

            // Store binary file
            std::filesystem::path output_path{arguments.out_path};
            output_path /= file_name.stem();
            output_path += ".minimiser";
            std::ofstream outfile{output_path, std::ios::binary};
            for (auto && hash : minimiser_table)
            {
                if (hash.second > cutoff)
                {
                    outfile.write(reinterpret_cast<const char*>(&hash.first), sizeof(hash.first));
                    ++count;
                }
            }

            // Store header file
            output_path = arguments.out_path;
            output_path /= file_name.stem();
            output_path += ".header";
            std::ofstream headerfile{output_path};
            headerfile << static_cast<uint64_t>(arguments.kmer_size) << '\t'
                        << arguments.window_size << '\t'
                        << cutoff << '\t'
                        << count << '\n';

            count = 0;
            minimiser_table.clear();
        }
    };

    call_parallel_on_bins(worker, arguments);
}

} // namespace raptor

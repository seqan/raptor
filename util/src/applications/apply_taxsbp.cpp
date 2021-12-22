// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <filesystem>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/algorithm/detail/execution_handler_parallel.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/utility/views/chunk.hpp>

#include <seqan3/core/debug_stream.hpp>

struct config
{
    uint64_t threads{1u};
    bool skip_gzip{false};

    std::filesystem::path input_directory{};
    std::filesystem::path output_directory{};
    std::filesystem::path taxsbp_binning_path{};
    std::filesystem::path genome_updater_accession_path{};
    std::filesystem::path assembly_summary_path{};
};

class positive_integer_validator
{
public:
    using option_value_type = size_t;

    positive_integer_validator() = default;
    positive_integer_validator(bool const is_zero_positive_) : is_zero_positive{is_zero_positive_} {}

    void operator() (option_value_type const & val) const
    {
        if (!is_zero_positive && !val)
        {
            throw seqan3::validation_error{"The value must be a positive integer."};
        }
    }

    std::string get_help_page_message () const
    {
        if (is_zero_positive)
            return "Value must be a positive integer or 0.";
        else
            return "Value must be a positive integer.";
    }

private:
    bool is_zero_positive{false};
};

struct char_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = char;
    using sequence_legal_alphabet = char;
};

inline std::unordered_map<std::string, std::filesystem::path> parse_assembly_summary(config const & cfg)
{
    using namespace std::literals;

    std::unordered_map<std::string, std::filesystem::path> assembly_accession_to_path{};

    std::ifstream file{cfg.assembly_summary_path};
    // std::string const extension = "_genomic.fna"s + (!cfg.skip_gzip ? ".gz"s : ""s);
    std::string const extension{"_genomic.fna.gz"}; // Actually gzipped by default

    if (file.is_open())
    {
        std::string line;
        std::string assembly_accession;
        std::filesystem::path assembly_filename;
        std::filesystem::path assembly_path;

        while (std::getline(file, line))
        {
            auto split_line = line
                            | std::views::split('\t')
                            | std::views::transform([] (auto && rng)
                              {
                                  return std::string_view(std::addressof(*rng.begin()), std::ranges::distance(rng));
                              });

            auto it = split_line.begin();
            assembly_accession = *it; // [0]
            std::advance(it, 19);
            assembly_filename = *it; // [19]
            assembly_filename = assembly_filename.filename();
            assembly_filename += extension;

            assembly_path = cfg.input_directory;
            assembly_path /= assembly_filename;

            assembly_accession_to_path[assembly_accession] = assembly_path;
        }
    }
    else
    {
        throw std::runtime_error{std::string{"Could not open file "} + cfg.assembly_summary_path.string()};
    }

    return assembly_accession_to_path;
}

inline std::unordered_map<std::string, std::string> parse_genome_updater_accession(config const & cfg)
{
    std::unordered_map<std::string, std::string> refseq_to_assembly_accession{};

    std::ifstream file{cfg.genome_updater_accession_path};

    if (file.is_open())
    {
        std::string line;
        std::string refseq_accession;
        std::string assembly_accession;

        while (std::getline(file, line))
        {
            auto split_line = line
                            | std::views::split('\t')
                            | std::views::transform([] (auto && rng)
                              {
                                  return std::string_view(std::addressof(*rng.begin()), std::ranges::distance(rng));
                              });

            auto it = split_line.begin();
            std::advance(it, 1);
            assembly_accession = *it; // [1]
            std::advance(it, 2);
            refseq_accession = *it; // [3]

            refseq_to_assembly_accession[refseq_accession] = assembly_accession;
        }
    }
    else
    {
        throw std::runtime_error{std::string{"Could not open file "} + cfg.genome_updater_accession_path.string()};
    }

    return refseq_to_assembly_accession;
}

inline std::unordered_map<uint64_t, std::vector<std::string>> parse_sbp_binning(config const & cfg)
{
    std::unordered_map<uint64_t, std::vector<std::string>> bin_to_refseq_accession{};

    std::ifstream file{cfg.taxsbp_binning_path};

    if (file.is_open())
    {
        std::string line;
        std::string refseq_accession;
        uint64_t bin_index;

        while (std::getline(file, line))
        {
            auto split_line = line
                            | std::views::split('\t')
                            | std::views::transform([] (auto && rng)
                              {
                                  return std::string_view(std::addressof(*rng.begin()), std::ranges::distance(rng));
                              });

            auto it = split_line.begin();
            refseq_accession = *it; // [0]
            std::advance(it, 5);
            std::from_chars((*it).begin(), (*it).end(), bin_index); // [4]

            bin_to_refseq_accession[bin_index].push_back(refseq_accession);
        }
    }
    else
    {
        throw std::runtime_error{std::string{"Could not open file "} + cfg.taxsbp_binning_path.string()};
    }

    return bin_to_refseq_accession;
}

inline void apply_taxsbp(config const & cfg)
{
    using namespace std::literals;

    std::unordered_map<std::string, std::filesystem::path> const assembly_accession_to_path{parse_assembly_summary(cfg)};
    std::unordered_map<std::string, std::string> const refseq_to_assembly_accession{parse_genome_updater_accession(cfg)};
    std::unordered_map<uint64_t, std::vector<std::string>> const bin_to_refseq_accession{parse_sbp_binning(cfg)};
    size_t const num_bins = bin_to_refseq_accession.size();
    size_t const n_zero = std::to_string(num_bins).length();

    auto worker = [&] (auto && chunk_view, auto &&)
        {
            for (uint64_t bin_index : chunk_view)
            {
                std::string bin_index_as_string = std::to_string(bin_index);
                std::string padded_bin_index = std::string(n_zero - bin_index_as_string.length(), '0') + bin_index_as_string;
                std::string filename = "bin_"s + padded_bin_index + ".fasta"s + (!cfg.skip_gzip ? ".gz"s : ""s);

                std::filesystem::path output_path = cfg.output_directory;
                output_path /= filename;

                seqan3::sequence_file_output output_file{output_path};
                output_file.options.fasta_blank_before_id = false;

                if (auto bin_to_refseq = bin_to_refseq_accession.find(bin_index); bin_to_refseq != bin_to_refseq_accession.end())
                {
                    for (auto && refseq_accession : bin_to_refseq->second)
                    {
                        if (auto refseq_to_assembly = refseq_to_assembly_accession.find(refseq_accession); refseq_to_assembly != refseq_to_assembly_accession.end())
                        {
                            if (auto assembly_to_path = assembly_accession_to_path.find(refseq_to_assembly->second); assembly_to_path != assembly_accession_to_path.end())
                            {
                                seqan3::sequence_file_input<char_traits, seqan3::fields<seqan3::field::seq, seqan3::field::id>> input_file{assembly_to_path->second};
                                for (auto && [seq, id] : input_file)
                                {
                                    if (id.substr(0, refseq_accession.size()) == refseq_accession)
                                    {
                                        output_file.emplace_back(seq, id);
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        };

    size_t const chunk_size = std::bit_ceil(num_bins / cfg.threads);
    auto chunked_view = std::views::iota(0u, num_bins) | seqan3::views::chunk(chunk_size);
    seqan3::detail::execution_handler_parallel executioner{cfg.threads};
    executioner.bulk_execute(std::move(worker), std::move(chunked_view), [](){});
}

int main(int argc, char ** argv)
{
    seqan3::argument_parser parser{"apply_taxsbp", argc, argv, seqan3::update_notifications::off};
    parser.info.author = "Enrico Seiler";
    parser.info.author = "enrico.seiler@fu-berlin.de";
    parser.info.short_description = "Split RefSeq according to taxsbp";
    parser.info.version = "0.0.1";

    config cfg{};

    parser.add_option(cfg.input_directory,
                      '\0',
                      "input",
                      "Provide the path to directory containing the RefSeq files.",
                      seqan3::option_spec::required,
                      seqan3::input_directory_validator{});

    parser.add_option(cfg.output_directory,
                      '\0',
                      "output",
                      "Provide an output directory.",
                      seqan3::option_spec::required,
                      seqan3::output_directory_validator{});

    parser.add_option(cfg.taxsbp_binning_path,
                      '\0',
                      "taxsbp",
                      "Provide the taxsbp binning file.",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{});

    parser.add_option(cfg.genome_updater_accession_path,
                      '\0',
                      "genome_update",
                      "Provide the genome_updater's updated_sequence_accession.txt file.",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{});

    parser.add_option(cfg.assembly_summary_path,
                      '\0',
                      "assembly_summary",
                      "Provide the assembly_summary file.",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{});

    parser.add_flag(cfg.skip_gzip,
                    '\0',
                    "skip_gzip",
                    "Whether to skip gzipping the output files.");

    parser.add_option(cfg.threads,
                      '\0',
                      "threads",
                      "Number of threads to use.",
                      seqan3::option_spec::standard,
                      positive_integer_validator{});

    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        std::cerr << "[Error] " << ext.what() << '\n';
        std::exit(-1);
    }

    std::filesystem::create_directory(cfg.output_directory);
    apply_taxsbp(cfg);
}

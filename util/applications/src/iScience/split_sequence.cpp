// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <filesystem>

#include <sharg/all.hpp>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>

#include <hibf/contrib/std/chunk_view.hpp>

static inline std::vector<std::string> sequence_extensions{
    seqan3::detail::valid_file_extensions<typename seqan3::sequence_file_input<>::valid_formats>()};

static inline std::vector<std::string> compression_extensions{[]()
                                                              {
                                                                  std::vector<std::string> result;
#ifdef SEQAN3_HAS_BZIP2
                                                                  result.push_back("bz2");
#endif
#ifdef SEQAN3_HAS_ZLIB
                                                                  result.push_back("gz");
                                                                  result.push_back("bgzf");
#endif
                                                                  return result;
                                                              }()}; // GCOVR_EXCL_LINE

static inline std::vector<std::string> combined_extensions{
    []()
    {
        if (compression_extensions.empty())
            return sequence_extensions; // GCOVR_EXCL_LINE
        std::vector<std::string> result;
        for (auto && sequence_extension : sequence_extensions)
        {
            result.push_back(sequence_extension);
            for (auto && compression_extension : compression_extensions)
                result.push_back(sequence_extension + std::string{'.'} + compression_extension);
        }
        return result;
    }()};

struct dna4_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

struct config
{
    uint64_t parts{};
    uint64_t length{};
    uint64_t sequence_size{};

    std::filesystem::path input_path{};
    std::filesystem::path output_path{};
};

class positive_integer_validator
{
public:
    using option_value_type = size_t;

    positive_integer_validator() = default;
    positive_integer_validator(bool const is_zero_positive_) : is_zero_positive{is_zero_positive_}
    {}

    void operator()(option_value_type const & val) const
    {
        if (!is_zero_positive && !val)
        {
            throw sharg::validation_error{"The value must be a positive integer."};
        }
    }

    std::string get_help_page_message() const
    {
        if (is_zero_positive)
            return "Value must be a positive integer or 0.";
        else
            return "Value must be a positive integer.";
    }

private:
    bool is_zero_positive{false};
};

inline void split_sequence(config const & cfg)
{
    seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>> input_sequence{cfg.input_path};

    std::string id;
    size_t part{};
    size_t const n_zero = std::to_string(cfg.parts).length();

    for (auto && split_sequence : (*input_sequence.begin()).sequence() | seqan::stl::views::chunk(cfg.length))
    {
        std::string part_as_string = std::to_string(part);
        std::string padded_parts = std::string(n_zero - part_as_string.length(), '0') + part_as_string;
        std::string filename = "bin_" + padded_parts + ".fa";
        std::filesystem::path out_path = cfg.output_path;
        out_path /= filename;
        seqan3::sequence_file_output output_sequence{out_path};
        output_sequence.options.fasta_blank_before_id = false;

        id = "bin_" + padded_parts;
        size_t const start = part++ * cfg.length;
        size_t const end = part * cfg.length;
        id += "_[" + std::to_string(start) + ',' + std::to_string(end) + ')';

        output_sequence.emplace_back(split_sequence, id);
    }
}

int main(int argc, char ** argv)
{
    sharg::parser parser{"split_sequence", argc, argv, sharg::update_notifications::off};
    parser.info.author = "Enrico Seiler";
    parser.info.author = "enrico.seiler@fu-berlin.de";
    parser.info.short_description = "Split a fasta into parts.";
    parser.info.version = "0.0.1";

    config cfg{};

    parser.add_option(cfg.input_path,
                      sharg::config{.short_id = '\0',
                                    .long_id = "input",
                                    .description = "Provide the path to a sequence to be split.",
                                    .required = true,
                                    .validator = sharg::input_file_validator{combined_extensions}});
    parser.add_option(cfg.output_path,
                      sharg::config{.short_id = '\0',
                                    .long_id = "output",
                                    .description = "Provide an output filepath. Default: directory of input."});
    parser.add_option(
        cfg.parts,
        sharg::config{.short_id = '\0',
                      .long_id = "parts",
                      .description =
                          "The number of parts to split the sequence into. Mutually exclusive with --length.",
                      .validator = positive_integer_validator{}});
    parser.add_option(cfg.length,
                      sharg::config{.short_id = '\0',
                                    .long_id = "length",
                                    .description = "The length of one part. Mutually exclusive with --parts.",
                                    .validator = positive_integer_validator{}});

    try
    {
        parser.parse();
    }
    catch (sharg::parser_error const & ext)
    {
        std::cerr << "[Error] " << ext.what() << '\n';
        std::exit(-1);
    }

    if (!parser.is_option_set("length") && !parser.is_option_set("parts"))
    {
        throw sharg::validation_error{"Set --length or --parts"};
    }

    // {
    //     seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>> input_sequence{cfg.input_path};
    //     auto sequence_count = std::ranges::distance(input_sequence);
    //     if (sequence_count == 0)
    //         throw sharg::validation_error{"Empty input sequence."};
    //     if (sequence_count > 1)
    //         throw sharg::validation_error{"Only one sequence per input sequence allowed."};
    // }

    if (!parser.is_option_set("length"))
    {
        seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>> input_sequence{cfg.input_path};
        cfg.sequence_size = std::ranges::distance((*input_sequence.begin()).sequence());
        cfg.length = cfg.sequence_size / cfg.parts;
    }

    if (!parser.is_option_set("parts"))
    {
        seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>> input_sequence{cfg.input_path};
        cfg.sequence_size = std::ranges::distance((*input_sequence.begin()).sequence());
        cfg.parts = cfg.sequence_size / cfg.length;
    }

    if (!parser.is_option_set("output"))
    {
        cfg.output_path = cfg.input_path.parent_path();
        if (cfg.output_path == "")
            cfg.output_path = ".";
    }

    sharg::output_directory_validator{}(cfg.output_path);

    split_sequence(cfg);
}

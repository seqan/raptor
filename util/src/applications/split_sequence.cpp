#include <seqan3/std/filesystem>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/utility/views/chunk.hpp>

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

inline void split_sequence(config const & cfg)
{
    seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>> input_sequence{cfg.input_path};

    std::string id;
    size_t part{};
    size_t const n_zero = std::to_string(cfg.parts).length();

    for (auto && split_sequence : (*input_sequence.begin()).sequence() | seqan3::views::chunk(cfg.length))
    {
        std::string part_as_string = std::to_string(part);
        std::string padded_parts = std::string(n_zero - part_as_string.length(), '0') + part_as_string;
        std::string filename = "bin_" + padded_parts + ".fasta";
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
    seqan3::argument_parser parser{"split_sequence", argc, argv, seqan3::update_notifications::off};
    parser.info.author = "Enrico Seiler";
    parser.info.author = "enrico.seiler@fu-berlin.de";
    parser.info.short_description = "Split a fasta into parts.";
    parser.info.version = "0.0.1";

    config cfg{};

    parser.add_option(cfg.input_path,
                      '\0',
                      "input",
                      "Provide the path to a sequence to be split.",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator<seqan3::sequence_file_input<>>{});

    parser.add_option(cfg.output_path,
                      '\0',
                      "output",
                      "Provide an output filepath. Default: directory of input.",
                      seqan3::option_spec::standard);

    parser.add_option(cfg.parts,
                      '\0',
                      "parts",
                      "The number of parts to split the sequence into. Mutually exclusive with --length.",
                      seqan3::option_spec::standard,
                      positive_integer_validator{});

    parser.add_option(cfg.length,
                      '\0',
                      "length",
                      "The length of one part. Mutually exclusive with --parts.",
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

    if (!parser.is_option_set("length") && !parser.is_option_set("parts"))
    {
        throw seqan3::validation_error{"Set --length or --parts"};
    }

    // {
    //     seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>> input_sequence{cfg.input_path};
    //     auto sequence_count = std::ranges::distance(input_sequence);
    //     if (sequence_count == 0)
    //         throw seqan3::validation_error{"Empty input sequence."};
    //     if (sequence_count > 1)
    //         throw seqan3::validation_error{"Only one sequence per input sequence allowed."};
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

    seqan3::output_directory_validator{}(cfg.output_path);

    split_sequence(cfg);
}

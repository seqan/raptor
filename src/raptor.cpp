#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/bit_manipulation.hpp>
#include <seqan3/range/views/async_input_buffer.hpp>

#include <raptor_build.hpp>
#include <raptor_search.hpp>
#include <shared.hpp>

void try_parsing(seqan3::argument_parser & parser)
{
    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        std::cerr << "[Error] " << ext.what() << '\n';
        std::exit(-1);
    }
}

struct power_of_two_validator
{
    using option_value_type = size_t;

    void operator() (option_value_type const & val) const
    {
        if (!seqan3::detail::is_power_of_two(val))
        {
            throw seqan3::validation_error{"The value must be a power of two."};
        }
    }

    std::string get_help_page_message () const
    {
        return "Value must be a power of two.";
    }
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

class size_validator
{
public:
    using option_value_type = std::string;

    size_validator(std::string const & pattern_) :
        pattern{pattern_}
    {}

    void operator()(option_value_type const & cmp) const
    {
        std::regex rgx(pattern);
        if (!std::regex_match(cmp, rgx))
            throw seqan3::validation_error{seqan3::detail::to_string("Value ", cmp, " must be an integer followed by [k,m,g,t] (case insensitive).")};
    }

    template <std::ranges::forward_range range_type>
        requires std::convertible_to<std::ranges::range_value_t<range_type>, option_value_type const &>
    void operator()(range_type const & v) const
    {
         std::for_each(v.begin(), v.end(), [&] (auto cmp) { (*this)(cmp); });
    }

    std::string get_help_page_message() const
    {
        return "Must be an integer followed by [k,m,g,t] (case insensitive).";
    }

private:
    std::string pattern;
};

class bin_validator
{
public:
    using option_value_type = std::vector<std::filesystem::path>;

    void operator() (option_value_type const & values) const
    {
        if (values.empty())
            throw seqan3::validation_error{"The list of input files cannot be empty."};

        for (auto && value : values)
        {
            try
            {
                sequence_file_validator(value);
            }
            catch (seqan3::validation_error const & exception)
            {
                if (value.extension() == ".minimiser")
                    minimiser_file_validator(value);
                else
                    throw exception;
            }
        }

        bool const is_minimiser_input = values[0].extension() == ".minimiser";

        for (auto && value : values)
            if (is_minimiser_input != (value.extension() == ".minimiser"))
                throw seqan3::validation_error{"You cannot mix sequence and minimiser files as input."};
    }

    std::string get_help_page_message() const
    {
        return "The input file must exist and read permissions must be granted. Valid file extensions are: [minimiser],"
               " or [embl, fasta, fa, fna, ffn, faa, frn, fas, fastq, fq, genbank, gb, gbk, sam] possibly followed by: "
               "[gz, bgzf, bz2].";
    }

private:
    seqan3::input_file_validator<seqan3::sequence_file_input<>> sequence_file_validator{};
    seqan3::input_file_validator<> minimiser_file_validator{{"minimiser"}};
};

inline void init_shared_meta(seqan3::argument_parser & parser)
{
    parser.info.app_name = "Raptor";
    parser.info.author = "Enrico Seiler";
    parser.info.citation = "Seiler, E. et al. (2020). Raptor: A fast and space-efficient pre-filter for"
                           " querying very large collections of nucleotide sequences. Manuscript in submission.";
    parser.info.date = "21-09-2020";
    parser.info.email = "enrico.seiler@fu-berlin.de";
    parser.info.long_copyright = R"(BSD 3-Clause License

Copyright (c) 2020, Enrico Seiler
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
    parser.info.short_description = "A fast and space-efficient pre-filter for querying very large collections of nucleotide sequences.";
    parser.info.url = "https://github.com/eseiler/raptor";
    parser.info.version = "1.0.0";
}

void init_top_level_parser(seqan3::argument_parser & parser)
{
    init_shared_meta(parser);
    parser.info.description.emplace_back("Binning Directories are a datastruture that can be used in various ways. "
                                         "What's a bin, how can it be used, etc.");

    parser.info.examples = {"./raptor build --help", "./raptor search --help"};
};

template <typename arguments_t>
    requires std::same_as<arguments_t, build_arguments> || std::same_as<arguments_t, search_arguments>
inline void init_shared_options(seqan3::argument_parser & parser, arguments_t & arguments)
{
    parser.add_option(arguments.window_size,
                      '\0',
                      "window",
                      "Choose the window size.",
                      seqan3::option_spec::DEFAULT,
                      positive_integer_validator{});
    parser.add_option(arguments.kmer_size,
                      '\0',
                      "kmer",
                      "Choose the kmer size.",
                      seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{1, 32});
    parser.add_option(arguments.threads,
                      '\0',
                      "threads",
                      "Choose the number of threads.",
                      seqan3::option_spec::DEFAULT,
                      positive_integer_validator{});
    parser.add_option(arguments.parts,
                      '\0',
                      "parts",
                      "Splits the index in this many parts.",
                      seqan3::option_spec::DEFAULT,
                      power_of_two_validator{});
}

inline void init_build_parser(seqan3::argument_parser & parser, build_arguments & arguments)
{
    init_shared_meta(parser);
    init_shared_options(parser, arguments);
    parser.add_positional_option(arguments.bin_path,
                                 "Provide a list of input files. One file per bin. ",
                                 bin_validator{});
    parser.add_option(arguments.out_path,
                      '\0',
                      "output",
                      "Provide an output filepath or an output directory if --compute-minimiser is used.",
                      seqan3::option_spec::REQUIRED);
    parser.add_option(arguments.size,
                      '\0',
                      "size",
                      "Choose the size of the resulting IBF.",
                      seqan3::option_spec::REQUIRED,
                      size_validator{"\\d+\\s{0,1}[k,m,g,t,K,M,G,T]"});
    parser.add_option(arguments.hash,
                      '\0',
                      "hash",
                      "Choose the number of hashes.",
                      seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{1, 5});
    parser.add_flag(arguments.compressed,
                    '\0',
                    "compressed",
                    "Build a compressed IBF.");
    parser.add_flag(arguments.compute_minimiser,
                    '\0',
                    "compute-minimiser",
                    "Computes minimisers using cutoffs from Mantis (Pandey et al.). Does not create the index.");
}

inline void init_search_parser(seqan3::argument_parser & parser, search_arguments & arguments)
{
    init_shared_meta(parser);
    init_shared_options(parser, arguments);
    parser.add_option(arguments.ibf_file,
                      '\0',
                      "index",
                      "Provide a valid path to an IBF. Parts: Without suffix _0",
                      seqan3::option_spec::REQUIRED,
                      seqan3::input_file_validator{});
    parser.add_option(arguments.query_file,
                      '\0',
                      "query",
                      "Provide a path to the query file.",
                      seqan3::option_spec::REQUIRED,
                      seqan3::input_file_validator<seqan3::sequence_file_input<>>{});
    parser.add_option(arguments.out_file,
                      '\0',
                      "output",
                      "Please provide a valid path to the output.",
                      seqan3::option_spec::REQUIRED,
                      seqan3::output_file_validator{});
    parser.add_option(arguments.errors,
                      '\0',
                      "error",
                      "Choose the number of errors.",
                      seqan3::option_spec::DEFAULT,
                      positive_integer_validator{true});
    parser.add_option(arguments.tau,
                      '\0',
                      "tau",
                      "Threshold for probabilistic models.",
                      seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{0, 1});
    parser.add_option(arguments.threshold,
                      '\0',
                      "threshold",
                      "If set, this threshold is used instead of the probabilistic models.",
                      seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{0, 1});
    parser.add_option(arguments.pattern_size,
                      '\0',
                      "pattern",
                      "Choose the pattern size. Default: Use median of sequence lengths in query file.");
    parser.add_flag(arguments.compressed,
                    '\0',
                    "compressed",
                    "Build a compressed IBF.");
    parser.add_flag(arguments.write_time,
                    '\0',
                    "time",
                    "Write timing file.",
                    seqan3::option_spec::ADVANCED);
}

void run_build(seqan3::argument_parser & parser)
{
    build_arguments arguments{};
    init_build_parser(parser, arguments);
    try_parsing(parser);

    // ==========================================
    // Various checks.
    // ==========================================

    arguments.bins = arguments.bin_path.size();

    if (parser.is_option_set("window"))
    {
        if (arguments.kmer_size > arguments.window_size)
            throw seqan3::argument_parser_error{"The k-mer size cannot be bigger than the window size."};
    }
    else
        arguments.window_size = arguments.kmer_size;

    if (parser.is_option_set("compute-minimiser"))
    {
        try
        {
            seqan3::output_directory_validator{}(arguments.out_path);
        }
        catch (seqan3::argument_parser_error const & ext)
        {
            std::cerr << "[Error] " << ext.what() << '\n';
            std::exit(-1);
        }
    }
    else
    {
        try
        {
            seqan3::output_file_validator{}(arguments.out_path);
        }
        catch (seqan3::argument_parser_error const & ext)
        {
            std::cerr << "[Error] " << ext.what() << '\n';
            std::exit(-1);
        }
    }

    // ==========================================
    // Process --size.
    // ==========================================
    arguments.size.erase(std::remove(arguments.size.begin(), arguments.size.end(), ' '), arguments.size.end());

    size_t multiplier{};

    switch (std::tolower(arguments.size.back()))
    {
        case 't':
            multiplier = 8ull * 1024ull * 1024ull * 1024ull * 1024ull;
            break;
        case 'g':
            multiplier = 8ull * 1024ull * 1024ull * 1024ull;
            break;
        case 'm':
            multiplier = 8ull * 1024ull * 1024ull;
            break;
        case 'k':
            multiplier = 8ull * 1024ull;
            break;
        default:
            throw seqan3::argument_parser_error{"Use {k, m, g, t} to pass size. E.g., --size 8g."};
    }

    size_t size{};
    std::from_chars(arguments.size.data(), arguments.size.data() + arguments.size.size() - 1, size);
    size *= multiplier;
    arguments.bits = size / (((arguments.bins + 63) >> 6) << 6);

    // ==========================================
    // Dispatch
    // ==========================================
    raptor_build(arguments);
};

void run_search(seqan3::argument_parser & parser)
{
    search_arguments arguments{};
    init_search_parser(parser, arguments);
    try_parsing(parser);

    // ==========================================
    // Various checks.
    // ==========================================

    if (parser.is_option_set("window"))
    {
        if (arguments.kmer_size > arguments.window_size)
            throw seqan3::argument_parser_error{"The k-mer size cannot be bigger than the window size."};
    }
    else
        arguments.window_size = arguments.kmer_size;

    arguments.treshold_was_set = parser.is_option_set("treshold");

    // ==========================================
    // Process --pattern.
    // ==========================================
    if (!arguments.pattern_size)
    {
        std::vector<uint64_t> sequence_lengths{};
        seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>> query_in{arguments.query_file};
        for (auto & [seq] : query_in | seqan3::views::async_input_buffer(16))
        {
            sequence_lengths.push_back(std::ranges::size(seq));
        }
        std::sort(sequence_lengths.begin(), sequence_lengths.end());
        arguments.pattern_size = sequence_lengths[sequence_lengths.size()/2];
    }

    // ==========================================
    // Dispatch
    // ==========================================
    raptor_search(arguments);
};

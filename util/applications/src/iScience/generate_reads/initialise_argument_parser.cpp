// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <sharg/all.hpp>

#include <raptor/argument_parsing/validators.hpp>

#include "cmd_arguments.hpp"

namespace sharg::custom
{

template <>
struct parsing<raptor::util::generate_reads::weight>
{
    static inline std::unordered_map<std::string_view, raptor::util::generate_reads::weight> const enumeration_names{
        {"from_file_sizes", raptor::util::generate_reads::weight::from_file_sizes},
        {"from_hll_sketches", raptor::util::generate_reads::weight::from_hll_sketches},
        {"from_kmer_counts", raptor::util::generate_reads::weight::from_kmer_counts},
        {"from_weight_column", raptor::util::generate_reads::weight::from_weight_column},
        {"from_uniform_distribution", raptor::util::generate_reads::weight::from_uniform_distribution}};
};

} // namespace sharg::custom

namespace raptor::util::generate_reads
{

class read_output_validator
{
private:
    static inline std::vector<std::string> fasta_fastq_extensions = []()
    {
        auto extensions = seqan3::format_fasta::file_extensions;
        std::ranges::copy(seqan3::format_fastq::file_extensions, std::back_inserter(extensions));
        return extensions;
    }();

    // Duplicated because static initialization order is not guaranteed.
    // raptor::detail::compression_extensions might not be initialized yet.
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
                return fasta_fastq_extensions; // GCOVR_EXCL_LINE
            std::vector<std::string> result;
            for (auto const & sequence_extension : fasta_fastq_extensions)
            {
                result.push_back(sequence_extension);
                for (auto const & compression_extension : compression_extensions)
                    result.push_back(sequence_extension + std::string{'.'} + compression_extension);
            }
            return result;
        }()};

public:
    using option_value_type = std::string;

    read_output_validator() = default;
    read_output_validator(read_output_validator const &) = default;
    read_output_validator & operator=(read_output_validator const &) = default;
    read_output_validator(read_output_validator &&) = default;
    read_output_validator & operator=(read_output_validator &&) = default;
    ~read_output_validator() = default;

    void operator()(option_value_type const & value) const
    {
        std::filesystem::path const out_path{value};
        std::filesystem::path const out_dir{out_path.parent_path()};
        if (!out_dir.empty())
        {
            std::error_code ec{};
            std::filesystem::create_directories(out_dir, ec);
            if (ec)
                throw sharg::validation_error{
                    sharg::detail::to_string("Failed to create directory \"", out_dir.c_str(), "\": ", ec.message())};
        }

        validator(out_path);
    }

    std::string get_help_page_message() const
    {
        return sharg::detail::to_string("Write permissions must be granted. Valid file extensions are ",
                                        fasta_fastq_extensions,
#if defined(SEQAN3_HAS_BZIP2) || defined(SEQAN3_HAS_ZLIB)
                                        ", possibly followed by ",
                                        raptor::detail::compression_extensions,
#endif
                                        ". ");
    }

private:
    sharg::output_file_validator validator{sharg::output_file_open_options::open_or_create, combined_extensions};
};

void initialise_argument_parser(sharg::parser & parser, cmd_arguments & arguments)
{
    parser.info.author = "Enrico Seiler";
    parser.info.author = "enrico.seiler@fu-berlin.de";
    parser.info.short_description = "Simulate reads from bins";
    parser.info.description.emplace_back("This utility can be used to simulate reads from bins.");
    parser.info.description.emplace_back("The number of reads simulated is controlled by weights. "
                                         "There are several ways to infer weights:");
    parser.info.description.emplace_back("1. From file sizes: The number of reads per bin is "
                                         "proportional to the size of the respectice bin file.");
    parser.info.description.emplace_back("2. [Default] From HLL sketches: The number of reads per bin is "
                                         "proportional to the HyperLogLog estimate for the respective bin. Canonical "
                                         "32-mers are used to create the sketch.");
    parser.info.description.emplace_back("3. From k-mer counts: The number of reads per bin is "
                                         "proportional to the number of canonical 32-mers in the respective bin. This "
                                         "approach may have a high memory consumption.");
    parser.info.description.emplace_back("4. From a weight column in the input file: The number of reads per bin is "
                                         "proportional to the weights specified in the input file.");
    parser.info.description.emplace_back("5. From a uniform distribution: The number of reads is equally distributed "
                                         "across all bins.");
    parser.info.version = "0.0.1";
    parser.info.examples = {"generate_reads --input weights.txt --weights from_uniform_distribution --output "
                            "reads.fasta --errors 2 --read_length 250 --number_of_reads 100000 --threads 4 "
                            "--print_weights"};

    parser.add_option(
        arguments.bin_file_path,
        sharg::config{
            .short_id = '\0',
            .long_id = "input",
            .description =
                "A file containing file names, one per line. If `--weights from_weight_column` is used, each line must "
                "additionally contain a weight (integer) after the file name. File name and weight must be separated "
                "by a tab. If any other weight mode is used, the weight column is ignored and may not even be present.",
            .required = true,
            .validator = sharg::input_file_validator{}});
    parser.add_list_item("", "Example file without weights:");
    parser.add_list_item("", "```");
    parser.add_list_item("", "/absolute/path/to/file1.fasta");
    parser.add_list_item("", "/absolute/path/to/file2.fa.gz");
    parser.add_list_item("", "```");
    parser.add_list_item("", "Example file with weights, note that column separator is a tab:");
    parser.add_list_item("", "```");
    parser.add_list_item("", "/absolute/path/to/file1.fasta\t27");
    parser.add_list_item("", "/absolute/path/to/file2.fa.gz\t78");
    parser.add_list_item("", "```");
    parser.add_option(arguments.weight_mode,
                      sharg::config{.short_id = '\0',
                                    .long_id = "weights",
                                    .description = "How to infer weights.",
                                    .validator = sharg::value_list_validator{
                                        std::views::values(sharg::custom::parsing<weight>::enumeration_names)}});
    parser.add_option(arguments.output_filename,
                      sharg::config{.short_id = '\0',
                                    .long_id = "output",
                                    .description = "Where to write the simulated reads",
                                    .required = true,
                                    .validator = read_output_validator{}});
    parser.add_list_item("",
                         "Using FASTA extensions will create a FASTA file (reads only), while FASTQ extensions will "
                         "create a FASTQ file (reads and quality). The quality is always set to the maximum value.");
    parser.add_option(arguments.errors,
                      sharg::config{.short_id = '\0',
                                    .long_id = "errors",
                                    .description = "The number of errors.",
                                    .validator = raptor::positive_integer_validator{true}});
    parser.add_option(arguments.read_length,
                      sharg::config{.short_id = '\0',
                                    .long_id = "read_length",
                                    .description = "The read length.",
                                    .validator = raptor::positive_integer_validator{}});
    parser.add_option(arguments.number_of_reads,
                      sharg::config{.short_id = '\0',
                                    .long_id = "number_of_reads",
                                    .description = "The number of reads.",
                                    .validator = raptor::positive_integer_validator{}});
    parser.add_option(arguments.threads,
                      sharg::config{.short_id = '\0',
                                    .long_id = "threads",
                                    .description = "Number of threads to use.",
                                    .validator = raptor::positive_integer_validator{}});
    parser.add_flag(arguments.print_weights,
                    sharg::config{.short_id = '\0',
                                  .long_id = "print_weights",
                                  .description = "Print weights and read counts for each bin."});
}

} // namespace raptor::util::generate_reads

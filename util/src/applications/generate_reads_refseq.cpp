#include <random>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/algorithm/detail/execution_handler_parallel.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/utility/views/chunk.hpp>

struct cmd_arguments
{
    std::filesystem::path bin_file_path{};
    std::vector<std::filesystem::path> bin_path{};
    std::filesystem::path output_directory{};
    uint8_t errors{2u};
    uint32_t read_length{100u};
    uint32_t number_of_reads{1ULL<<20};
    uint64_t threads{1u};
};

struct char_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = char;
    using sequence_legal_alphabet = char;
};

struct dna4_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
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

inline size_t count_records_in_fasta(std::filesystem::path const & filename)
{
    seqan3::sequence_file_input<char_traits, seqan3::fields<seqan3::field::id>> fin{filename};
    return std::ranges::distance(fin);
}

void run_program(cmd_arguments const & arguments)
{
    std::uniform_int_distribution<uint32_t> read_error_position_dis(0, arguments.read_length - 1);
    std::uniform_int_distribution<uint8_t> dna4_rank_dis(0, 3);

    size_t const number_of_bins{arguments.bin_path.size()};
    uint32_t const reads_per_bin = arguments.number_of_reads / number_of_bins;

    std::vector<seqan3::phred42> const quality(arguments.read_length, seqan3::assign_rank_to(40u, seqan3::phred42{}));

    auto worker = [&] (auto && zipped_view, auto &&)
    {
        for (auto && [bin_file, bin_number] : zipped_view)
        {
            std::mt19937_64 rng(bin_number);
            uint32_t read_counter{bin_number * reads_per_bin};
            // Immediately invoked initialising lambda expession (IIILE).
            std::filesystem::path const out_file = [&] ()
                                                   {
                                                       std::filesystem::path out_file = arguments.output_directory;
                                                       if (bin_file.extension() == ".gz")
                                                           out_file /= bin_file.stem().stem();
                                                       else
                                                           out_file /= bin_file.stem();
                                                       out_file += ".fastq";
                                                       return out_file;
                                                    }();

            size_t const number_of_records{count_records_in_fasta(bin_file)};
            uint32_t const reads_per_record = (reads_per_bin + number_of_records - 1) / number_of_records;

            seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>> fin{bin_file};
            seqan3::sequence_file_output fout{out_file};

            uint16_t bin_read_counter{};

            for (auto const & [seq] : fin)
            {
                uint64_t const reference_length = std::ranges::size(seq);
                std::uniform_int_distribution<uint64_t> read_start_dis(0, reference_length - arguments.read_length);
                for (uint32_t current_read_number = 0;
                    current_read_number < reads_per_record && bin_read_counter < reads_per_bin;
                    ++current_read_number, ++read_counter, ++bin_read_counter)
                {
                    uint64_t const read_start_pos = read_start_dis(rng);
                    std::vector<seqan3::dna4> read = seq
                                                   | seqan3::views::slice(read_start_pos, read_start_pos + arguments.read_length)
                                                   | seqan3::views::to<std::vector>;

                    for (uint8_t error_count = 0; error_count < arguments.errors; ++error_count)
                    {
                        uint32_t const error_pos = read_error_position_dis(rng);
                        seqan3::dna4 const current_base = read[error_pos];
                        seqan3::dna4 new_base = current_base;
                        while (new_base == current_base)
                            seqan3::assign_rank_to(dna4_rank_dis(rng), new_base);
                        read[error_pos] = new_base;
                    }

                    fout.emplace_back(read, std::to_string(read_counter), quality);
                }
            }
        }
    };

    size_t const chunk_size = std::bit_ceil(number_of_bins / arguments.threads);
    auto chunked_view = seqan3::views::zip(arguments.bin_path, std::views::iota(0u)) | seqan3::views::chunk(chunk_size);
    seqan3::detail::execution_handler_parallel executioner{arguments.threads};
    executioner.bulk_execute(std::move(worker), std::move(chunked_view), [](){});
}

void initialise_argument_parser(seqan3::argument_parser & parser, cmd_arguments & arguments)
{
    parser.info.author = "Enrico Seiler";
    parser.info.author = "enrico.seiler@fu-berlin.de";
    parser.info.short_description = "Generate reads from bins.";
    parser.info.version = "0.0.1";
    parser.info.examples = {"./generate_reads_refseq --output ./reads_e2 all_bins.txt"};
    parser.add_positional_option(arguments.bin_file_path,
                                 "Provide a path to a file containing one path to a bin per line.",
                                 seqan3::input_file_validator{});
    parser.add_option(arguments.output_directory,
                      '\0',
                      "output",
                      "Provide the base dir where the reads should be written to.",
                      seqan3::option_spec::required,
                      seqan3::output_directory_validator{});
    parser.add_option(arguments.errors,
                      '\0',
                      "errors",
                      "The number of errors.",
                      seqan3::option_spec::standard,
                      positive_integer_validator{true});
    parser.add_option(arguments.read_length,
                      '\0',
                      "read_length",
                      "The read length.",
                      seqan3::option_spec::standard,
                      positive_integer_validator{});
    parser.add_option(arguments.number_of_reads,
                      '\0',
                      "number_of_reads",
                      "The number of reads.",
                      seqan3::option_spec::standard,
                      positive_integer_validator{});
    parser.add_option(arguments.threads,
                      '\0',
                      "threads",
                      "Number of threads to use.",
                      seqan3::option_spec::standard,
                      positive_integer_validator{});
}

int main(int argc, char ** argv)
{
    seqan3::argument_parser myparser{"build_ibf", argc, argv, seqan3::update_notifications::off};
    cmd_arguments arguments{};
    initialise_argument_parser(myparser, arguments);
    try
    {
         myparser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        std::cout << "[Error] " << ext.what() << "\n";
        return -1;
    }

    std::ifstream istrm{arguments.bin_file_path};
    std::string line;
    seqan3::input_file_validator validator{};

    while (std::getline(istrm, line))
    {
        if (!line.empty())
        {
            std::filesystem::path bin_path{line};
            validator(bin_path);
            arguments.bin_path.push_back(std::move(bin_path));
        }
    }

    size_t const number_of_bins{arguments.bin_path.size()};

    if (arguments.errors > arguments.read_length)
        throw seqan3::argument_parser_error{"Cannot have more errors than the read is long."};

    if (number_of_bins > arguments.number_of_reads)
        throw seqan3::argument_parser_error{"Must simulate at least one read per bin."};

    if (arguments.number_of_reads % number_of_bins)
        throw seqan3::argument_parser_error{"The number of reads must distribute evenly over the bins."};

    std::filesystem::create_directory(arguments.output_directory);
    run_program(arguments);

    return 0;
}

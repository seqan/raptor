#include <seqan3/std/filesystem>
#include <mutex>
#include <unordered_set>
#include <vector>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/algorithm/detail/execution_handler_parallel.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/utility/views/chunk.hpp>

inline constexpr static uint64_t adjust_seed(uint8_t const kmer_size, uint64_t const seed = 0x8F3F73B5CF1C9ADEULL) noexcept
{
    return seed >> (64u - 2u * kmer_size);
}

struct dna4_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

struct config
{
    uint32_t window_size{};
    uint8_t kmer_size{};
    uint8_t threads{1u};

    std::vector<std::filesystem::path> bin_path{};
    std::filesystem::path out_path{};
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

inline void compute_minimisers(config const & cfg)
{
    auto minimiser_view = seqan3::views::minimiser_hash(seqan3::ungapped{cfg.kmer_size},
                                                        seqan3::window_size{cfg.window_size},
                                                        seqan3::seed{adjust_seed(cfg.kmer_size)});

    std::vector<uint64_t> minimiser_counts;
    minimiser_counts.reserve(cfg.bin_path.size());
    std::mutex push_back_mutex;
    std::unordered_set<uint64_t> all_minimiser_set{};
    std::mutex merge_mutex;

    auto worker = [&] (auto && file_range, auto &&)
        {
            std::unordered_set<uint64_t> minimiser_set{};

            for (auto && file_name : file_range)
            {
                seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>> fin{file_name};

                for (auto & [seq] : fin)
                    for (auto && hash : seq | minimiser_view)
                        minimiser_set.insert(hash);

                {
                    std::lock_guard const lock{push_back_mutex};
                    minimiser_counts.emplace_back(minimiser_set.size());
                }
                {
                    std::lock_guard const lock{merge_mutex};
                    all_minimiser_set.merge(minimiser_set);
                }

                minimiser_set.clear();
            }
        };

    size_t const chunk_size = std::ceil<size_t>(cfg.bin_path.size() / cfg.threads);
    auto chunked_view = cfg.bin_path | seqan3::views::chunk(chunk_size);
    seqan3::detail::execution_handler_parallel executioner{cfg.threads};
    executioner.bulk_execute(worker, std::move(chunked_view), [](){});


    std::ofstream output{cfg.out_path.string()};
    output << all_minimiser_set.size() << '\n';
    for (auto const & count : minimiser_counts)
        output << count << '\n';
}

int main(int argc, char ** argv)
{
    seqan3::argument_parser parser{"count_minimiser", argc, argv, seqan3::update_notifications::off};
    parser.info.author = "Enrico Seiler";
    parser.info.author = "enrico.seiler@fu-berlin.de";
    parser.info.short_description = "Count minimiser.";
    parser.info.version = "0.0.1";

    config cfg{};

    parser.add_positional_option(cfg.bin_path,
                                 "Provide a list of input files (one file per bin).",
                                 seqan3::input_file_validator<seqan3::sequence_file_input<>>{});

    parser.add_option(cfg.out_path,
                      '\0',
                      "output",
                      "Provide an output filepath.",
                      seqan3::option_spec::required,
                      seqan3::output_file_validator{seqan3::output_file_open_options::create_new});

    parser.add_option(cfg.window_size,
                      '\0',
                      "window",
                      "Choose the window size.",
                      seqan3::option_spec::required,
                      positive_integer_validator{});

    parser.add_option(cfg.kmer_size,
                      '\0',
                      "kmer",
                      "Choose the kmer size.",
                      seqan3::option_spec::required,
                      seqan3::arithmetic_range_validator{1, 32});

    parser.add_option(cfg.threads,
                      '\0',
                      "threads",
                      "Choose the number of threads.",
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

    compute_minimisers(cfg);
}

#include <seqan3/std/charconv>
#include <limits> // Needed for gcc11 robin_hood

#include <robin_hood.h>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/algorithm/detail/execution_handler_parallel.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/utility/views/chunk.hpp>

#include <shared.hpp>

template <std::copy_constructible algorithm_t>
void call_parallel_on_bins(algorithm_t && worker, build_arguments const & arguments)
{
    size_t const chunk_size = std::clamp<size_t>(std::bit_ceil(arguments.bins / arguments.threads),
                                                 8u,
                                                 64u);
    auto chunked_view = seqan3::views::zip(arguments.bin_path, std::views::iota(0u)) |
                        seqan3::views::chunk(chunk_size);
    seqan3::detail::execution_handler_parallel executioner{arguments.threads};
    executioner.bulk_execute(std::move(worker), std::move(chunked_view), [](){});
}

template <bool compressed>
class ibf_factory
{
public:
    ibf_factory(build_arguments const & args) : arguments{args} {}

    ibf_factory() = default;
    ibf_factory(ibf_factory const &) = default;
    ibf_factory(ibf_factory &&) = default;
    ibf_factory & operator=(ibf_factory const &) = default;
    ibf_factory & operator=(ibf_factory &&) = default;
    ~ibf_factory() = default;

    template <typename view_t = std::ranges::empty_view<int>>
        requires (!compressed)
    auto ibf(view_t && hash_filter_view = std::ranges::empty_view<int>()) const
    {
        return construct(std::move(hash_filter_view));
    }

    template <typename view_t = std::ranges::empty_view<int>>
        requires compressed
    auto ibf(view_t && hash_filter_view = std::ranges::empty_view<int>()) const
    {
        auto tmp = construct(std::move(hash_filter_view));

        return seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed>{std::move(tmp)};
    }

private:
    build_arguments const arguments{};

    template <typename view_t = std::ranges::empty_view<int>>
    auto construct(view_t && hash_filter_view = std::ranges::empty_view<int>()) const
    {
        using sequence_file_t = seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>>;

        seqan3::interleaved_bloom_filter<> ibf{seqan3::bin_count{arguments.bins},
                                               seqan3::bin_size{arguments.bits / arguments.parts},
                                               seqan3::hash_function_count{arguments.hash}};

        if constexpr (std::same_as<view_t, std::ranges::empty_view<int>>)
        {
            auto worker = [&] (auto && zipped_view, auto &&)
            {
                auto hash_view = seqan3::views::minimiser_hash(seqan3::ungapped{arguments.kmer_size},
                                                               seqan3::window_size{arguments.window_size},
                                                               seqan3::seed{adjust_seed(arguments.kmer_size)});

                for (auto && [file_name, bin_number] : zipped_view)
                    for (auto && [seq] : sequence_file_t{file_name})
                        for (auto && value : seq | hash_view)
                            ibf.emplace(value, seqan3::bin_index{bin_number});
            };

            call_parallel_on_bins(worker, arguments);
        }
        else
        {
            auto worker = [&] (auto && zipped_view, auto &&)
            {
                auto hash_view = seqan3::views::minimiser_hash(seqan3::ungapped{arguments.kmer_size},
                                                               seqan3::window_size{arguments.window_size},
                                                               seqan3::seed{adjust_seed(arguments.kmer_size)})
                                 | hash_filter_view;

                for (auto && [file_name, bin_number] : zipped_view)
                    for (auto && [seq] : sequence_file_t{file_name})
                        for (auto && value : seq | hash_view)
                            ibf.emplace(value, seqan3::bin_index{bin_number});
            };

            call_parallel_on_bins(worker, arguments);
        }

        return ibf;
    }
};

template <bool compressed>
void run_program(build_arguments const & arguments)
{
    ibf_factory<compressed> generator{arguments};

    if (arguments.parts == 1u)
    {
        auto ibf = generator.ibf();
        std::ofstream os{arguments.out_path, std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(ibf);
    }
    else
    {
        std::vector<std::vector<size_t>> association(arguments.parts);
        size_t next_power_of_four{4u};

        if (arguments.parts == 4u) // one-to-one
        {
            for (size_t i : std::views::iota(0u, arguments.parts))
                association[i] = std::vector<size_t>{i};
        }
        else if (arguments.parts == 2u) // More than 1 prefix per part
        {
            association[0] = std::vector<size_t>{0, 1};
            association[1] = std::vector<size_t>{2, 3};
        }
        else // More parts than prefixes
        {
            // How long must the suffix be such that 4^suffix_length >= arguments.parts
            size_t suffix_length{0};
            for (; 0b100 << (2 * suffix_length) < arguments.parts; ++suffix_length) {}
            next_power_of_four = 0b100 << (2 * suffix_length);

            size_t const prefixes_per_part = next_power_of_four / arguments.parts;

            for (size_t i : std::views::iota(0u, next_power_of_four))
                association[i/prefixes_per_part].push_back(i);
        }

        for (size_t part : std::views::iota(0u, arguments.parts))
        {
            size_t const mask{next_power_of_four - 1};
            auto filter_view = std::views::filter([&] (auto && hash)
                { return std::ranges::find(association[part], hash & mask) != association[part].end(); });

            auto ibf = generator.ibf(filter_view);
            std::filesystem::path out_path{arguments.out_path};
            out_path += "_" + std::to_string(part);
            std::ofstream os{out_path, std::ios::binary};
            cereal::BinaryOutputArchive oarchive{os};
            oarchive(ibf);
        }
    }
}

inline bool check_for_fasta_format(std::vector<std::string> const & valid_extensions, std::string const & file_path)
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

inline void compute_minimisers(build_arguments const & arguments)
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
            uint16_t cutoff{default_cutoff};

            for (auto && [file_name, bin_number] : zipped_view)
            {
                (void) bin_number;
                seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>> fin{file_name};

                for (auto & [seq] : fin)
                    for (auto && hash : seq | minimiser_view)
                        minimiser_table[hash] = std::min<uint8_t>(254u, minimiser_table[hash] + 1);
                        // The hash table stores how often a minimiser appears. It does not matter whether a minimiser appears
                        // 50 times or 2000 times, it is stored regardless because the biggest cutoff value is 50. Hence,
                        // the hash table stores only values up to 254 to save memory.

                // Since the curoffs are based on the filesize of a gzipped fastq file, we try account for the other cases:
                // We multiply by two if we have fasta input.
                // We divide by 3 if the input is not compressed.
                bool const is_compressed = file_name.extension() == ".gz" || file_name.extension() == ".bgzf" || file_name.extension() == ".bz2";
                bool const is_fasta = is_compressed ? check_for_fasta_format(seqan3::format_fasta::file_extensions, file_name.stem())
                                                    : check_for_fasta_format(seqan3::format_fasta::file_extensions, file_name.extension());
                size_t const filesize = std::filesystem::file_size(file_name) * (is_fasta ? 2 : 1) / (is_compressed ? 1 : 3);

                for (size_t k = 0; k < cutoff_bounds.size(); ++k)
                {
                    if (filesize <= cutoff_bounds[k])
                    {
                        cutoff = cutoffs[k];
                        break;
                    }
                }

                // Store binary file
                std::ofstream outfile{arguments.out_path.string() + file_name.stem().string() + ".minimiser", std::ios::binary};
                for (auto && hash : minimiser_table)
                {
                    if (hash.second > cutoff)
                    {
                        outfile.write(reinterpret_cast<const char*>(&hash.first), sizeof(hash.first));
                        ++count;
                    }
                }

                // Store header file
                std::ofstream headerfile{arguments.out_path.string() + file_name.stem().string() + ".header"};
                headerfile << static_cast<uint64_t>(arguments.kmer_size) << '\t'
                           << arguments.window_size << '\t'
                           << cutoff << '\t'
                           << count << '\n';

                count = 0;
                cutoff = default_cutoff;
                minimiser_table.clear();
            }
        };

    call_parallel_on_bins(worker, arguments);
}

template <bool compressed>
void build_from_binary(build_arguments const & arguments)
{
    seqan3::interleaved_bloom_filter<> ibf{seqan3::bin_count{arguments.bins},
                                            seqan3::bin_size{arguments.bits / arguments.parts},
                                            seqan3::hash_function_count{arguments.hash}};

    auto worker = [&] (auto && zipped_view, auto &&)
        {
            uint64_t read_number;

            for (auto && [file_name, bin_number] : zipped_view)
            {
                std::ifstream infile{file_name, std::ios::binary};

                while(infile.read(reinterpret_cast<char*>(&read_number), sizeof(read_number)))
                    ibf.emplace(read_number, seqan3::bin_index{bin_number});
            }
        };

    call_parallel_on_bins(std::move(worker), arguments);

    if constexpr (compressed)
    {
        seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> cibf{std::move(ibf)};

        std::ofstream os{arguments.out_path, std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(cibf);
    }
    else
    {
        std::ofstream os{arguments.out_path, std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(ibf);
    }
}

void raptor_build(build_arguments const & arguments)
{
    if (arguments.compute_minimiser)
    {
        compute_minimisers(arguments);
        return;
    }

    if (arguments.bin_path[0].extension() == ".minimiser")
    {
        if (arguments.compressed)
            build_from_binary<true>(arguments);
        else
            build_from_binary<false>(arguments);
        return;
    }

    if (arguments.compressed)
        run_program<true>(arguments);
    else
        run_program<false>(arguments);
    return;
}

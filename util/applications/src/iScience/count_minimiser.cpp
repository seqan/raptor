// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <filesystem>
#include <mutex>
#include <vector>

#include <sharg/all.hpp>

#include <seqan3/core/algorithm/detail/execution_handler_parallel.hpp>

#include <raptor/argument_parsing/parse_bin_path.hpp>
#include <raptor/argument_parsing/validators.hpp>
#include <raptor/file_reader.hpp>

#include <hibf/contrib/robin_hood.hpp>
#include <hibf/contrib/std/chunk_view.hpp>
#include <hibf/contrib/std/join_with_view.hpp>
#include <hibf/contrib/std/zip_view.hpp>

struct config
{
    uint32_t window_size{};
    uint8_t kmer_size{};
    uint8_t threads{1u};

    std::filesystem::path bin_file{};
    std::vector<std::vector<std::string>> bin_path{};
    std::filesystem::path out_path{};
};

void compute_minimisers(config const & cfg)
{
    std::vector<uint64_t> minimiser_counts(cfg.bin_path.size());
    robin_hood::unordered_node_set<uint64_t> all_minimiser_set{};
    std::mutex merge_mutex;

    auto worker = [&](auto && zipped_view, auto &&)
    {
        robin_hood::unordered_flat_set<uint64_t> minimiser_set{};
        raptor::file_reader<raptor::file_types::sequence> reader{seqan3::ungapped{cfg.kmer_size}, cfg.window_size};

        for (auto && [file_names, bin_number] : zipped_view)
        {
            reader.hash_into(file_names, std::inserter(minimiser_set, minimiser_set.begin()));

            minimiser_counts[bin_number] = minimiser_set.size();

            {
                std::lock_guard const lock{merge_mutex};
                all_minimiser_set.insert(minimiser_set.begin(), minimiser_set.end());
            }

            minimiser_set.clear();
        }
    };

    size_t const chunk_size = std::max<size_t>(1, std::floor(cfg.bin_path.size() / static_cast<double>(cfg.threads)));
    auto chunked_view =
        seqan::stl::views::zip(cfg.bin_path, std::views::iota(0u)) | seqan::stl::views::chunk(chunk_size);
    seqan3::detail::execution_handler_parallel executioner{cfg.threads};
    executioner.bulk_execute(worker, std::move(chunked_view), []() {});

    std::ofstream output{cfg.out_path.string()};
    output << "# k-mer size: " << static_cast<size_t>(cfg.kmer_size) << '\n';
    output << "# Window size: " << cfg.window_size << '\n';
    output << "# Unique k-mers across all: " << all_minimiser_set.size() << '\n';
    output << "# File\tUnique k-mer count\n";
    for (auto && [file_names, count] : seqan::stl::views::zip(cfg.bin_path, minimiser_counts))
    {
        for (auto && elem : seqan::stl::views::join_with(file_names, ' '))
            output << elem;
        output << '\t';
        output << count << '\n';
    }
}

int main(int argc, char ** argv)
{
    sharg::parser parser{"count_minimiser", argc, argv, sharg::update_notifications::off};
    parser.info.author = "Enrico Seiler";
    parser.info.author = "enrico.seiler@fu-berlin.de";
    parser.info.short_description = "Count minimiser.";
    parser.info.version = "0.0.1";

    config cfg{};

    parser.add_positional_option(
        cfg.bin_file,
        sharg::config{.description = "File containing file names.", .validator = sharg::input_file_validator{}});

    parser.add_option(
        cfg.out_path,
        sharg::config{.short_id = '\0',
                      .long_id = "output",
                      .description = "Provide an output filepath.",
                      .required = true,
                      .validator = sharg::output_file_validator{sharg::output_file_open_options::create_new}});
    parser.add_option(cfg.window_size,
                      sharg::config{.short_id = '\0',
                                    .long_id = "window",
                                    .description = "Choose the window size.",
                                    .required = true,
                                    .validator = raptor::positive_integer_validator{}});
    parser.add_option(cfg.kmer_size,
                      sharg::config{.short_id = '\0',
                                    .long_id = "kmer",
                                    .description = "Choose the kmer size.",
                                    .required = true,
                                    .validator = sharg::arithmetic_range_validator{1, 32}});
    parser.add_option(cfg.threads,
                      sharg::config{.short_id = '\0',
                                    .long_id = "threads",
                                    .description = "Choose the number of threads.",
                                    .validator = raptor::positive_integer_validator{}});

    try
    {
        parser.parse();
    }
    catch (sharg::parser_error const & ext)
    {
        std::cerr << "[Error] " << ext.what() << '\n';
        std::exit(-1);
    }

    raptor::detail::parse_bin_path(cfg.bin_file, cfg.bin_path, false);

    compute_minimisers(cfg);
}

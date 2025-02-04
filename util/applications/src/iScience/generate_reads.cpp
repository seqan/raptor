// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <sharg/all.hpp>

#include "generate_reads/cmd_arguments.hpp"

namespace raptor::util::generate_reads
{

void apply_weights(std::vector<size_t> & vec, size_t const number_of_reads);
void print_weights(std::vector<size_t> const & vec, std::string_view const first_line);
void infer_weights_from_file_size(cmd_arguments & arguments);
void infer_weights_from_hll_sketches(cmd_arguments & arguments);
void infer_weights_from_kmer_counts(cmd_arguments & arguments);
void infer_weights_from_uniform_distribution(cmd_arguments & arguments);
void initialise_argument_parser(sharg::parser & parser, cmd_arguments & arguments);
void parse_input(cmd_arguments & arguments);
size_t simulate_reads(cmd_arguments const & arguments);

} // namespace raptor::util::generate_reads

namespace raptor::util::generate_reads
{

void run(std::vector<std::string> const & cli_arguments)
{
    sharg::parser parser{"generate_reads", cli_arguments, sharg::update_notifications::off};
    cmd_arguments arguments{};
    initialise_argument_parser(parser, arguments);
    parser.parse();

    parse_input(arguments);

    switch (arguments.weight_mode)
    {
    case weight::from_hll_sketches:
        infer_weights_from_hll_sketches(arguments);
        break;
    case weight::from_file_sizes:
        infer_weights_from_file_size(arguments);
        break;
    case weight::from_kmer_counts:
        infer_weights_from_kmer_counts(arguments);
        break;
    case weight::from_weight_column:
        break;
    case weight::from_uniform_distribution:
        infer_weights_from_uniform_distribution(arguments);
        break;
    default:
        throw sharg::parser_error{"Unknown weight mode."};
    }

    if (arguments.print_weights)
        print_weights(arguments.number_of_reads_per_bin, "Weights:");

    if (arguments.errors > arguments.read_length)
        throw sharg::parser_error{"Cannot have more errors than the read is long."};

    if (arguments.number_of_bins > arguments.number_of_reads)
        throw sharg::parser_error{"Must simulate at least one read per bin."};

    if (arguments.number_of_bins != arguments.number_of_reads_per_bin.size())
        throw sharg::parser_error{"arguments.number_of_bins (" + std::to_string(arguments.number_of_bins)
                                  + " != arguments.number_of_reads_per_bin.size()"
                                  + std::to_string(arguments.number_of_reads_per_bin.size()) + ")"};

    apply_weights(arguments.number_of_reads_per_bin, arguments.number_of_reads);

    if (arguments.print_weights)
        print_weights(arguments.number_of_reads_per_bin, "Read counts:");

    size_t const total_reads = simulate_reads(arguments);

    std::cout << "\nTotal simulated reads: " << total_reads << '\n';
}

} // namespace raptor::util::generate_reads

int main(int argc, char ** argv)
{
    try
    {
        // Use the (user's environment) default locale for output. Usually results in numbers having separators.
        std::cout.imbue(std::locale(""));
        raptor::util::generate_reads::run({argv, argv + argc});
    }
    catch (sharg::parser_error const & ext)
    {
        std::cerr << "[Error] " << ext.what() << '\n';
        std::exit(-1);
    }

    return 0;
}

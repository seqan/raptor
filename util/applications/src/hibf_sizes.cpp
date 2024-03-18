// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <iostream>

#include <sharg/detail/to_string.hpp>
#include <sharg/exceptions.hpp>
#include <sharg/parser.hpp>

#include <raptor/index.hpp>
#include <raptor/search/load_index.hpp>

struct config
{
    std::filesystem::path input{};
    std::filesystem::path output{"sizes.txt"};
};

inline void set_up_parser(sharg::parser & parser, config & cfg)
{
    parser.info.version = "1.0.0";
    parser.info.author = "Svenja Mehringer";
    parser.info.email = "svenja.mehringer@fu-berlin.de";
    parser.info.short_description = "Read IBF sizes from an HIBF";

    parser.add_subsection("Main options:");
    parser.add_option(cfg.input,
                      sharg::config{.short_id = '\0',
                                    .long_id = "index",
                                    .description = "The input must be an HIBF index computed via raptor layout/build. ",
                                    .required = true,
                                    .validator = sharg::input_file_validator{}});
    parser.add_option(cfg.output,
                      sharg::config{.short_id = '\0',
                                    .long_id = "output",
                                    .description = "The output. ",
                                    .required = true,
                                    .validator = sharg::output_file_validator{}});
}

int main(int argc, char const * argv[])
{
    sharg::parser parser{"hibf_sizes", argc, argv, sharg::update_notifications::off};
    parser.info.version = "1.0.0";

    config cfg{};
    set_up_parser(parser, cfg);

    try
    {
        parser.parse();
    }
    catch (sharg::parser_error const & ext)
    {
        std::cerr << "[ERROR] " << ext.what() << '\n';
        return -1;
    }

    raptor::raptor_index<raptor::index_structure::hibf> index{};
    raptor::detail::load_index(index, cfg.input);
    auto const & next_ibf_id = index.ibf().next_ibf_id;
    auto const & ibf_vector = index.ibf().ibf_vector;
    assert(next_ibf_id.size() == ibf_vector.size());
    size_t const number_of_ibfs = next_ibf_id.size();

    std::vector<size_t> ibf_levels = [&]()
    {
        std::vector<size_t> current_ibf_levels(number_of_ibfs);
        std::vector<size_t> old_ibf_levels(number_of_ibfs, 1u);
        while (current_ibf_levels != old_ibf_levels)
        {
            old_ibf_levels = current_ibf_levels;
            for (size_t ibf_idx = 0; ibf_idx < number_of_ibfs; ++ibf_idx)
            {
                for (size_t const next_ibf_idx : next_ibf_id[ibf_idx])
                {
                    if (next_ibf_idx != ibf_idx) // there is a lower level for this merged tb
                        current_ibf_levels[next_ibf_idx] = current_ibf_levels[ibf_idx] + 1;
                }
            }
        }
        return current_ibf_levels;
    }();

    size_t const number_of_levels = std::ranges::max(ibf_levels) + 1u;
    std::vector<size_t> size_per_level(number_of_levels);
    for (size_t ibf_idx = 0; ibf_idx < number_of_ibfs; ++ibf_idx)
        size_per_level[ibf_levels[ibf_idx]] += ibf_vector[ibf_idx].bit_size();

    std::ofstream output{cfg.output};
    if (!output.good() || !output.is_open())
        throw std::logic_error{"Could not open file " + cfg.output.string() + " for reading"};

    output << "IDX\tLEVEL\tSIZE\n";
    for (size_t ibf_idx = 0; ibf_idx < number_of_ibfs; ++ibf_idx)
    {
        output << ibf_idx << '\t' << ibf_levels[ibf_idx] << '\t' << ibf_vector[ibf_idx].bit_size() << '\n';
    }

    std::cout << "LEVEL\tSIZE\n";
    for (size_t i = 0; i < number_of_levels; ++i)
    {
        std::cout << i << '\t' << size_per_level[i] << '\n';
    }
}

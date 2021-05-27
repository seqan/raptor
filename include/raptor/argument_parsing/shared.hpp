#pragma once

#include <raptor/argument_parsing/validators.hpp>
#include <raptor/shared.hpp>

namespace raptor
{

void init_shared_meta(seqan3::argument_parser & parser);
void try_parsing(seqan3::argument_parser & parser);

template <typename arguments_t>
void init_shared_options(seqan3::argument_parser & parser, arguments_t & arguments)
{
    static_assert(std::same_as<arguments_t, build_arguments> || std::same_as<arguments_t, search_arguments>);

    parser.add_option(arguments.threads,
                      '\0',
                      "threads",
                      "Choose the number of threads.",
                      seqan3::option_spec::standard,
                      positive_integer_validator{});
    parser.add_option(arguments.parts,
                      '\0',
                      "parts",
                      "Splits the index in this many parts.",
                      seqan3::option_spec::standard,
                      power_of_two_validator{});
}

} // namespace raptor

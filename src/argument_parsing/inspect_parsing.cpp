// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::inspect_parsing.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <raptor/argument_parsing/inspect_parsing.hpp>
#include <raptor/argument_parsing/validators.hpp>
#include <raptor/inspect/inspect.hpp>

namespace raptor
{

void init_inspect_parser(sharg::parser & parser, inspect_arguments & arguments, std::string & format)
{
    parser.info.short_description = "Prints the metadata of a Raptor index";
    parser.info.description.emplace_back(
        "Prints the parameters a Raptor index was built with, e.g., the window size, the shape, and the false "
        "positive rate. For an HIBF, also prints its layout parameters and how far updates have changed it.");
    parser.info.description.emplace_back(
        "By default, only the parameters and the file names stored at the beginning of the index are read. This is "
        "fast even for very large indices. Use --detailed to load the whole index and report the structure of its "
        "IBFs.");
    parser.info.examples.emplace_back("raptor inspect --index raptor.index");
    parser.info.examples.emplace_back("raptor inspect --index raptor.index --user-bins --output metadata.txt");
    parser.info.examples.emplace_back("raptor inspect --index raptor.index --detailed --format json");
    parser.info.synopsis.emplace_back(
        "raptor inspect --index <file> [--output <file>] [--format text|json] [--user-bins] [--detailed]");

    parser.add_subsection("General options");
    parser.add_option(arguments.index_file,
                      sharg::config{.short_id = '\0',
                                    .long_id = "index",
                                    .description = "Provide a valid path to an index. Parts: Without suffix _0",
                                    .required = true});
    parser.add_option(arguments.out_file,
                      sharg::config{.short_id = '\0',
                                    .long_id = "output",
                                    .description = "Write the metadata to this file instead of the standard output.",
                                    .default_message = "Standard output",
                                    .validator = output_file_validator{}});
    // A string rather than raptor::inspect_format: sharg lists the allowed values only for a value_list_validator.
    parser.add_option(format,
                      sharg::config{.short_id = '\0',
                                    .long_id = "format",
                                    .description = "The output format. \\fBtext\\fP prints the same \"## Key = "
                                                   "value\" lines as the header of a raptor search result. "
                                                   "\\fBjson\\fP prints a single JSON object.",
                                    .validator = sharg::value_list_validator{"text", "json"}});

    parser.add_subsection("Content options");
    parser.add_flag(arguments.user_bins,
                    sharg::config{.short_id = '\0',
                                  .long_id = "user-bins",
                                  .description = "List the files of every user bin, one line per user bin."});
    parser.add_flag(arguments.detailed,
                    sharg::config{.short_id = '\0',
                                  .long_id = "detailed",
                                  .description = "Load the whole index to also report the structure of its IBFs. "
                                                 "Takes as long and as much memory as loading the index for a "
                                                 "search; a partitioned index is loaded one part at a time."});
}

void inspect_parsing(sharg::parser & parser)
{
    inspect_arguments arguments{};
    std::string format{"text"};
    init_inspect_parser(parser, arguments, format);
    parser.parse();

    arguments.format = format == "json" ? inspect_format::json : inspect_format::text;

    std::filesystem::path const partitioned_index_file = arguments.index_file.string() + "_0";
    bool const index_is_monolithic = std::filesystem::exists(arguments.index_file);
    arguments.is_partitioned = std::filesystem::exists(partitioned_index_file);

    if (index_is_monolithic && arguments.is_partitioned)
    {
        throw sharg::validation_error{sharg::detail::to_string("Ambiguous index. Both monolithic (",
                                                               arguments.index_file.c_str(),
                                                               ") and partitioned index (",
                                                               partitioned_index_file.c_str(),
                                                               ") exist. Please rename the monolithic index.")};
    }

    // The remaining parts are checked once the number of parts has been read from the index.
    sharg::input_file_validator{}(arguments.is_partitioned ? partitioned_index_file : arguments.index_file);

    raptor_inspect(arguments);
}

} // namespace raptor

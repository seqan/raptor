// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::raptor_inspect.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <fstream>
#include <iostream>

#include <sharg/exceptions.hpp>

#include <raptor/inspect/inspect.hpp>

#include "index_metadata.hpp"
#include "write_metadata.hpp"

namespace raptor
{

void raptor_inspect(inspect_arguments const & arguments)
{
    // Read before opening the output file: A failure must not leave an empty or truncated output behind.
    detail::index_metadata const metadata = detail::read_index_metadata(arguments);

    std::ofstream file{};
    if (!arguments.out_file.empty())
        file.open(arguments.out_file);
    std::ostream & stream = arguments.out_file.empty() ? std::cout : file;

    switch (arguments.format)
    {
    case inspect_format::text:
        detail::write_text(stream, metadata, arguments.user_bins);
        break;
    case inspect_format::json:
        detail::write_json(stream, metadata, arguments.user_bins);
        break;
    }

    stream.flush();
    if (!stream)
        throw sharg::parser_error{"Failed to write the metadata."}; // GCOVR_EXCL_LINE
}

} // namespace raptor

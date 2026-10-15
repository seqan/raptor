// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::raptor_inspect.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <raptor/argument_parsing/inspect_arguments.hpp>

namespace raptor
{

/*!\brief Writes the metadata of an index.
 * \param[in] arguments The inspect arguments.
 * \details
 * Reads the parameters of the index at `arguments.index_file` without loading its data and writes them to
 * `arguments.out_file`, or to `std::cout` if it is empty. With `arguments.detailed`, loads the whole index to also
 * report its structure.
 */
void raptor_inspect(inspect_arguments const & arguments);

} // namespace raptor

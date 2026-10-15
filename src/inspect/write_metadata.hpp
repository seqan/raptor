// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::detail::write_text and raptor::detail::write_json.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <iosfwd>

#include "index_metadata.hpp"

namespace raptor::detail
{

/*!\brief Writes `metadata` as `## Key = value` lines, like the header of a `raptor search` result.
 * \param[out] stream The stream to write to.
 * \param[in] metadata The metadata to write.
 * \param[in] user_bins Whether to list the files of every user bin, as `#<id>\t<file>[,<file>...]`, directly after
 *                      the number of user bins.
 * \details
 * If `metadata.structure` is set, the output ends with one tab-separated line per IBF, preceded by the column
 * names in a line starting with `#IBF`.
 */
void write_text(std::ostream & stream, index_metadata const & metadata, bool const user_bins);

/*!\brief Writes `metadata` as a single JSON object.
 * \param[out] stream The stream to write to.
 * \param[in] metadata The metadata to write.
 * \param[in] user_bins Whether to list the files of every user bin.
 * \details Contains the same values as raptor::detail::write_text.
 */
void write_json(std::ostream & stream, index_metadata const & metadata, bool const user_bins);

} // namespace raptor::detail

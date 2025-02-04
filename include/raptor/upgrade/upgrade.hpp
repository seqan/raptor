// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::raptor_upgrade.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <raptor/argument_parsing/upgrade_arguments.hpp>

namespace raptor
{

void raptor_upgrade(upgrade_arguments & arguments);

} // namespace raptor

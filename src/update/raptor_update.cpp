// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::raptor_update.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <fstream> // for basic_ifstream, ifstream
#include <string>  // for basic_string
#include <utility> // for move
#include <vector>  // for vector

#include <cereal/archives/binary.hpp> // for BinaryInputArchive

#include <raptor/argument_parsing/update_arguments.hpp> // for update_arguments
#include <raptor/build/store_index.hpp>                 // for store_index
#include <raptor/index.hpp>                             // for raptor_index, hibf
#include <raptor/update/delete_user_bins.hpp>           // for delete_user_bins
#include <raptor/update/insert_user_bin.hpp>            // for insert_user_bin
#include <raptor/update/update.hpp>                     // for raptor_update

namespace raptor
{

void raptor_update(update_arguments const & arguments)
{
    std::ifstream index_file{arguments.index_file};
    cereal::BinaryInputArchive archive{index_file};
    raptor::raptor_index<index_structure::hibf> index;
    archive(index);

    // dump_index(index);
    if (!arguments.user_bins_to_delete.empty())
    {
        delete_user_bins(arguments, index);
        // dump_index(index);
    }
    if (!arguments.user_bins_to_insert.empty())
    {
        insert_user_bin(arguments, index);
        // dump_index(index);
    }

    store_index(arguments.out_path, std::move(index));
}

} // namespace raptor

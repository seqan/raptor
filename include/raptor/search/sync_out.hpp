// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::sync_out.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <filesystem>
#include <fstream>
#include <mutex>

#include <hibf/contrib/std/join_with_view.hpp>

#include <raptor/argument_parsing/search_arguments.hpp>

namespace raptor
{

class sync_out
{
public:
    sync_out() = default;
    sync_out(sync_out const &) = delete;             // std::ofstream
    sync_out & operator=(sync_out const &) = delete; // std::ofstream
    sync_out(sync_out &&) = delete;                  // std::mutex
    sync_out & operator=(sync_out &&) = delete;      // std::mutex
    ~sync_out() = default;

    sync_out(search_arguments const & arguments) : file{arguments.out_file}
    {}

    template <typename t>
    void write(t && data)
    {
        std::lock_guard<std::mutex> lock{write_mutex};
        file << std::forward<t>(data);
    }

    bool write_header(search_arguments const & arguments, size_t const hash_function_count)
    {
        file << "### Minimiser parameters\n";
        file << "## Window size = " << arguments.window_size << '\n';
        file << "## Shape = " << arguments.shape.to_string() << '\n';
        file << "## Shape size (length) = " << static_cast<uint16_t>(arguments.shape_size) << '\n';
        file << "## Shape count (number of 1s) = " << static_cast<uint16_t>(arguments.shape_weight) << '\n';
        file << "### Search parameters\n";
        file << "## Query file = " << arguments.query_file << '\n';
        file << "## Pattern size = " << arguments.query_length << '\n';
        file << "## Output file = " << arguments.out_file << '\n';
        file << "## Threads = " << static_cast<uint16_t>(arguments.threads) << '\n';
        file << "## tau = " << arguments.tau << '\n';
        file << "## p_max = " << arguments.p_max << '\n';
        file << "## Percentage threshold = " << arguments.threshold << '\n';
        file << "## Errors = " << static_cast<uint16_t>(arguments.errors) << '\n';
        file << "## Cache thresholds = " << std::boolalpha << arguments.cache_thresholds << '\n';
        file << "### Index parameters\n";
        file << "## Index = " << arguments.index_file << '\n';
        file << "## Index hashes = " << hash_function_count << '\n';
        file << "## Index parts = " << static_cast<uint16_t>(arguments.parts) << '\n';
        file << "## False positive rate = " << arguments.fpr << '\n';
        file << "## Index is HIBF = " << std::boolalpha << arguments.is_hibf << '\n';

        size_t user_bin_id{};
        for (auto const & file_list : arguments.bin_path)
        {
            file << '#' << user_bin_id << '\t';
            for (auto const elem : seqan::stl::views::join_with(file_list, ','))
                file << elem;
            file << '\n';
            ++user_bin_id;
        }

        file << "#QUERY_NAME\tUSER_BINS\n";

        return true;
    }

private:
    std::ofstream file;
    std::mutex write_mutex;
};

} // namespace raptor

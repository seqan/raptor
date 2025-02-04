// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <filesystem>

#include <sharg/exceptions.hpp>

#include <seqan3/core/debug_stream/detail/to_string.hpp>

inline void check_output_file(std::filesystem::path const & output_file)
{
    std::filesystem::path const output_directory = output_file.parent_path();
    std::error_code ec{};
    std::filesystem::create_directories(output_directory, ec);

    if (!output_directory.empty() && ec)
        throw sharg::parser_error{
            seqan3::detail::to_string("Failed to create directory\"", output_directory.c_str(), "\": ", ec.message())};
}

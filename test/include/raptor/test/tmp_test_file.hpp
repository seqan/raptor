// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <filesystem>  // for path
#include <ostream>     // for basic_ios, basic_ostream, ofstream
#include <string_view> // for basic_string_view, string_view

#include <seqan3/test/sandboxed_path.hpp> // for sandboxed_path, operator/
#include <seqan3/test/tmp_directory.hpp>  // for tmp_directory

namespace raptor::test
{

class tmp_test_file
{
public:
    tmp_test_file() = default;
    tmp_test_file(tmp_test_file const &) = delete;
    tmp_test_file & operator=(tmp_test_file const &) = delete;
    tmp_test_file(tmp_test_file &&) = delete;
    tmp_test_file & operator=(tmp_test_file &&) = delete;
    ~tmp_test_file() = default;

    template <typename... t>
    std::filesystem::path create(std::string_view const filename, t... content) const
    {
        std::filesystem::path const path = test_directory.path() / filename;
        std::ofstream file{path};
        (void)(file << ... << content); // (void) silences "state has no effect" warning for sizeof...(t) == 0
        return path;
    }

    std::filesystem::path path() const
    {
        return test_directory.path();
    }

private:
    seqan3::test::tmp_directory const test_directory{};
};

} // namespace raptor::test

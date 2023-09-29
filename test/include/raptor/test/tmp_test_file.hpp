// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <seqan3/test/tmp_directory.hpp>

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

// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <raptor/file_reader.hpp>
#include <raptor/test/expect_throw_msg.hpp>

struct noop_iterator
{
    using difference_type = std::ptrdiff_t;

    constexpr noop_iterator & operator*() noexcept
    {
        return *this;
    }

    constexpr noop_iterator & operator=(auto &&) noexcept
    {
        return *this;
    }

    constexpr noop_iterator & operator++() noexcept
    {
        return *this;
    }

    constexpr noop_iterator operator++(int) noexcept
    {
        return *this;
    }
};
static_assert(std::output_iterator<noop_iterator, uint64_t>);

TEST(file_reader, sequence_no_exist)
{
    raptor::file_reader<raptor::file_types::sequence> const reader{};
    std::vector<std::string> const filenames{"does_not_exist.fa"};
    std::string const & filename{filenames.front()};
    noop_iterator target{};
    auto predicate = [](uint64_t)
    {
        return true;
    };
    auto callback = [](uint64_t) {};

    constexpr std::string_view expected{"Could not open file does_not_exist.fa for reading."};
    EXPECT_THROW_MSG(reader.hash_into(filenames, target), seqan3::file_open_error, expected);
    EXPECT_THROW_MSG(reader.hash_into(filename, target), seqan3::file_open_error, expected);
    EXPECT_THROW_MSG(reader.hash_into_if(filenames, target, predicate), seqan3::file_open_error, expected);
    EXPECT_THROW_MSG(reader.hash_into_if(filename, target, predicate), seqan3::file_open_error, expected);
    EXPECT_THROW_MSG(reader.for_each_hash(filenames, callback), seqan3::file_open_error, expected);
    EXPECT_THROW_MSG(reader.for_each_hash(filename, callback), seqan3::file_open_error, expected);
}

TEST(file_reader, minimiser_no_exist)
{
    raptor::file_reader<raptor::file_types::minimiser> const reader{};
    std::vector<std::string> const filenames{"does_not_exist.minimiser"};
    std::string const & filename{filenames.front()};
    noop_iterator target{};
    auto predicate = [](uint64_t)
    {
        return true;
    };
    auto callback = [](uint64_t) {};

    constexpr std::string_view expected{"Failed to open file: \"does_not_exist.minimiser\": No such file or directory"};
    EXPECT_THROW_MSG(reader.hash_into(filenames, target), std::system_error, expected);
    EXPECT_THROW_MSG(reader.hash_into(filename, target), std::system_error, expected);
    EXPECT_THROW_MSG(reader.hash_into_if(filenames, target, predicate), std::system_error, expected);
    EXPECT_THROW_MSG(reader.hash_into_if(filename, target, predicate), std::system_error, expected);
    EXPECT_THROW_MSG(reader.for_each_hash(filenames, callback), std::system_error, expected);
    EXPECT_THROW_MSG(reader.for_each_hash(filename, callback), std::system_error, expected);
}

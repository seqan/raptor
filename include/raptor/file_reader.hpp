// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides raptor::file_reader.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

#include <raptor/adjust_seed.hpp>
#include <raptor/dna4_traits.hpp>

namespace raptor
{

enum class file_types
{
    sequence,
    minimiser
};

template <file_types file_type>
class file_reader
{};

template <>
class file_reader<file_types::sequence>
{
public:
    file_reader() = default;
    file_reader(file_reader const &) = default;
    file_reader(file_reader &&) = default; // GCOVR_EXCL_LINE
    file_reader & operator=(file_reader const &) = default;
    file_reader & operator=(file_reader &&) = default;
    ~file_reader() = default;

    explicit file_reader(seqan3::shape const shape, uint32_t const window_size) :
        minimiser_view{seqan3::views::minimiser_hash(shape,
                                                     seqan3::window_size{window_size},
                                                     seqan3::seed{adjust_seed(shape.count())})}
    {}

    template <std::output_iterator<uint64_t> it_t>
    void hash_into(std::vector<std::string> const & filenames, it_t target) const
    {
        for (auto && filename : filenames)
            hash_into(filename, target);
    }

    template <std::output_iterator<uint64_t> it_t>
    void hash_into(std::string const & filename, it_t target) const
    {
        sequence_file_t fin{filename};
        for (auto && record : fin)
            std::ranges::copy(record.sequence() | minimiser_view, target);
    }

    template <std::output_iterator<uint64_t> it_t>
    void hash_into_if(std::vector<std::string> const & filenames, it_t target, auto && pred) const
    {
        for (auto && filename : filenames)
            hash_into_if(filename, target, pred);
    }

    template <std::output_iterator<uint64_t> it_t>
    void hash_into_if(std::string const & filename, it_t target, auto && pred) const
    {
        sequence_file_t fin{filename};
        for (auto && record : fin)
            std::ranges::copy_if(record.sequence() | minimiser_view, target, pred);
    }

    void on_hash(std::vector<std::string> const & filenames, auto && callback) const
    {
        for (auto && filename : filenames)
            on_hash(filename, callback);
    }

    void on_hash(std::string const & filename, auto && callback) const
    {
        sequence_file_t fin{filename};
        for (auto && record : fin)
            callback(record.sequence() | minimiser_view);
    }

    void for_each_hash(std::vector<std::string> const & filenames, auto && callback) const
    {
        for (auto && filename : filenames)
            for_each_hash(filename, callback);
    }

    void for_each_hash(std::string const & filename, auto && callback) const
    {
        sequence_file_t fin{filename};
        for (auto && record : fin)
            std::ranges::for_each(record.sequence() | minimiser_view, callback);
    }

private:
    using sequence_file_t = seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>>;
    using view_t = decltype(seqan3::views::minimiser_hash(seqan3::shape{}, seqan3::window_size{}, seqan3::seed{}));
    view_t minimiser_view = seqan3::views::minimiser_hash(seqan3::shape{}, seqan3::window_size{}, seqan3::seed{});
};

template <>
class file_reader<file_types::minimiser>
{
public:
    file_reader() = default;
    file_reader(file_reader const &) = default;
    file_reader(file_reader &&) = default;
    file_reader & operator=(file_reader const &) = default;
    file_reader & operator=(file_reader &&) = default;
    ~file_reader() = default;

    explicit file_reader(seqan3::shape const, uint32_t const)
    {}

    template <std::output_iterator<uint64_t> it_t>
    void hash_into(std::vector<std::string> const & filenames, it_t target) const
    {
        for (auto && filename : filenames)
            hash_into(filename, target);
    }

    template <std::output_iterator<uint64_t> it_t>
    void hash_into(std::string const & filename, it_t target) const
    {
        std::ifstream fin{filename, std::ios::binary};
        uint64_t value;
        while (fin.read(reinterpret_cast<char *>(&value), sizeof(value)))
        {
            *target = value;
            ++target;
        }
    }

    template <std::output_iterator<uint64_t> it_t>
    void hash_into_if(std::vector<std::string> const & filenames, it_t target, auto && pred) const
    {
        for (auto && filename : filenames)
            hash_into_if(filename, target, pred);
    }

    template <std::output_iterator<uint64_t> it_t>
    void hash_into_if(std::string const & filename, it_t target, auto && pred) const
    {
        std::ifstream fin{filename, std::ios::binary};
        uint64_t value;
        while (fin.read(reinterpret_cast<char *>(&value), sizeof(value)))
            if (pred(value))
            {
                *target = value;
                ++target;
            }
    }

    void for_each_hash(std::vector<std::string> const & filenames, auto && callback) const
    {
        for (auto && filename : filenames)
            for_each_hash(filename, callback);
    }

    void for_each_hash(std::string const & filename, auto && callback) const
    {
        std::ifstream fin{filename, std::ios::binary};
        uint64_t value;
        while (fin.read(reinterpret_cast<char *>(&value), sizeof(value)))
            callback(value);
    }
};

} // namespace raptor

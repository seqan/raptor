// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides raptor::emplace_iterator.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <hibf/interleaved_bloom_filter.hpp>

namespace raptor
{

class emplace_iterator
{
public:
    using iterator_category = std::output_iterator_tag;
    using value_type = void;
    using difference_type = ptrdiff_t;
    using pointer = void;
    using reference = void;

    emplace_iterator() = delete;
    emplace_iterator(emplace_iterator const &) = default;
    emplace_iterator(emplace_iterator &&) = default;
    emplace_iterator & operator=(emplace_iterator const &) = default;
    emplace_iterator & operator=(emplace_iterator &&) = default;
    ~emplace_iterator() = default;

    explicit constexpr emplace_iterator(seqan::hibf::interleaved_bloom_filter & ibf, seqan::hibf::bin_index const idx) :
        ibf{std::addressof(ibf)},
        index{std::move(idx)}
    {}

    /* constexpr */ emplace_iterator & operator=(uint64_t const value) noexcept
    {
        assert(ibf != nullptr);
        ibf->emplace(std::move(value), index);
        return *this;
    }

    [[nodiscard]] constexpr emplace_iterator & operator*() noexcept
    {
        return *this;
    }

    constexpr emplace_iterator & operator++() noexcept
    {
        return *this;
    }

    constexpr emplace_iterator operator++(int) noexcept
    {
        return *this;
    }

private:
    seqan::hibf::interleaved_bloom_filter * ibf{nullptr};
    seqan::hibf::bin_index index{};
};

[[nodiscard]] inline constexpr emplace_iterator emplacer(seqan::hibf::interleaved_bloom_filter & ibf,
                                                         seqan::hibf::bin_index const idx)
{
    return emplace_iterator{ibf, std::move(idx)};
}

} // namespace raptor

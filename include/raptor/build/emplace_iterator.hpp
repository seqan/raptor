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

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

namespace raptor
{

template <typename container_t>
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

    explicit constexpr emplace_iterator(container_t & cont, seqan3::bin_index const idx) :
        container{std::addressof(cont)},
        index{std::move(idx)}
    {}

    constexpr emplace_iterator & operator=(uint64_t const value)
    {
        container->emplace(std::move(value), index);
        return *this;
    }

    [[nodiscard]] constexpr emplace_iterator & operator*()
    {
        return *this;
    }

    constexpr emplace_iterator & operator++()
    {
        return *this;
    }

    constexpr emplace_iterator operator++(int)
    {
        return *this;
    }

private:
    container_t * container{};
    seqan3::bin_index index{};
};

template <typename container_t>
[[nodiscard]] inline constexpr emplace_iterator<container_t> emplacer(container_t & cont, seqan3::bin_index const idx)
{
    return emplace_iterator<container_t>(cont, std::move(idx));
}

} // namespace raptor

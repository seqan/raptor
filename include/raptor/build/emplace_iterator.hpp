// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::emplace_iterator.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <cassert>  // for assert
#include <cstddef>  // for ptrdiff_t
#include <cstdint>  // for uint64_t
#include <iterator> // for output_iterator_tag
#include <memory>   // for addressof

#include <hibf/interleaved_bloom_filter.hpp> // for bin_index, interleaved_bloom_filter

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
        index{idx}
    {}

    /* constexpr */ emplace_iterator & operator=(uint64_t const value) noexcept
    {
        assert(ibf != nullptr);
        ibf->emplace(value, index);
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
    return emplace_iterator{ibf, idx};
}

} // namespace raptor

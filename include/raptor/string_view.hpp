// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <string_view>

//LCOV_EXCL_START
namespace raptor::detail
{

#if defined(__GNUC__) && (__GNUC__ < 10)
    template <typename iterator_t, typename sentinel_t>
    constexpr std::string_view string_view(iterator_t && begin, sentinel_t && end) noexcept
    {
        assert(end - begin >= 0);
        return std::string_view{&(*begin), static_cast<size_t>(end - begin)};
    }

    template <typename CharT, typename Traits>
    constexpr bool ends_with(std::basic_string_view<CharT, Traits> sv, std::basic_string_view<CharT, Traits> suffix) noexcept
    {
        return sv.size() >= suffix.size() && sv.compare(sv.size() - suffix.size(), std::basic_string_view<CharT, Traits>::npos, suffix) == 0;
    }
#else
    template <typename iterator_t, typename sentinel_t>
    constexpr std::string_view string_view(iterator_t && begin, sentinel_t && end) noexcept
    {
        return std::string_view{std::forward<iterator_t>(begin), std::forward<sentinel_t>(end)};
    }

    template <typename CharT, typename Traits>
    constexpr bool ends_with(std::basic_string_view<CharT, Traits> sv, std::basic_string_view<CharT, Traits> suffix) noexcept
    {
        return sv.ends_with(suffix);
    }
#endif

} // namespace raptor::detail
//LCOV_EXCL_STOP

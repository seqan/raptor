// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::cutoff.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <algorithm>   // for __equal, __find_if, equal, find_if
#include <array>       // for array
#include <cctype>      // for tolower
#include <cstddef>     // for size_t
#include <cstdint>     // for uint8_t, uint64_t
#include <filesystem>  // for path, operator==, file_size
#include <string>      // for basic_string, string
#include <string_view> // for basic_string_view, string_view
#include <vector>      // for vector

#include <seqan3/io/sequence_file/format_fasta.hpp> // for format_fasta

#include <hibf/misc/unreachable.hpp> // for unreachable

#include <raptor/argument_parsing/prepare_arguments.hpp> // for prepare_arguments

namespace raptor
{

class cutoff
{
public:
    cutoff() = default;
    cutoff(cutoff const &) = default;
    cutoff & operator=(cutoff const &) = default;
    cutoff(cutoff &&) = default;
    cutoff & operator=(cutoff &&) = default;
    ~cutoff() = default;

    cutoff(prepare_arguments const & arguments) : fixed_cutoff{arguments.kmer_count_cutoff}
    {
        if (arguments.use_filesize_dependent_cutoff)
        {
            cutoff_kind = cutoff_kinds::filesize_dependent;
        }
        else
        {
            cutoff_kind = cutoff_kinds::fixed;
        }
    }

    uint8_t get(std::filesystem::path const & filename) const noexcept
    {
        switch (cutoff_kind)
        {
        case cutoff_kinds::filesize_dependent:
            return impl(filename);
        case cutoff_kinds::fixed:
            return fixed_cutoff;
        default: // GCOVR_EXCL_LINE
            seqan::hibf::unreachable();
        }
    }

    static inline bool file_is_compressed(std::filesystem::path const & filepath) noexcept
    {
        std::filesystem::path const extension = filepath.extension();
        return extension == ".gz" || extension == ".bgzf" || extension == ".bz2";
    }

private:
    enum class cutoff_kinds : uint8_t
    {
        fixed,
        filesize_dependent
    };

    cutoff_kinds cutoff_kind{cutoff_kinds::fixed};

    uint8_t fixed_cutoff{};
    // Cutoffs and bounds from Mantis
    // Mantis ignores k-mers which appear less than a certain cutoff. The cutoff is based on the file size of a
    // gzipped fastq file. Small files have only a cutoff of 1 while big files have a cutoff value of 50.
    // https://doi.org/10.1016/j.cels.2018.05.021
    // Supplement Table S1
    // https://www.cell.com/cms/10.1016/j.cels.2018.05.021/attachment/0a3d402b-8b90-42c0-a709-22f246fd1759/mmc1.pdf
    static constexpr std::array<uint8_t, 4> const cutoffs{1u, 3u, 10u, 20u};
    static constexpr std::array<uint64_t, 4> const cutoff_bounds{314'572'800ULL,
                                                                 524'288'000ULL,
                                                                 1'073'741'824ULL,
                                                                 3'221'225'472ULL};

    uint8_t impl(std::filesystem::path const & filename) const noexcept
    {
        bool const is_compressed = file_is_compressed(filename);
        bool const is_fasta = check_for_fasta_format(filename);

        // Since the curoffs are based on the filesize of a gzipped fastq file, we try account for the other cases:
        // We multiply by two if we have fasta input.
        // We divide by 3 if the input is not compressed.
        size_t const filesize = std::filesystem::file_size(filename) * (is_fasta ? 2 : 1) / (is_compressed ? 1 : 3);

        uint8_t cutoff{50u};
        for (size_t i = 0; i < cutoff_bounds.size(); ++i)
        {
            if (filesize <= cutoff_bounds[i])
            {
                cutoff = cutoffs[i];
                break;
            }
        }

        return cutoff;
    }

    // NOLINTNEXTLINE(bugprone-exception-escape)
    static inline bool check_for_fasta_format(
        std::filesystem::path const & filepath,
        std::vector<std::string> const & valid_extensions = seqan3::format_fasta::file_extensions) noexcept
    {
        std::string const extension = file_is_compressed(filepath) ? filepath.stem() : filepath.extension();

        auto case_insensitive_string_ends_with = [&](std::string_view str, std::string_view suffix)
        {
            size_t const suffix_length{suffix.size()};
            size_t const str_length{str.size()};
            return suffix_length > str_length ? false
                                              : std::ranges::equal(str.substr(str_length - suffix_length),
                                                                   suffix,
                                                                   [](char const chr1, char const chr2)
                                                                   {
                                                                       return std::tolower(chr1) == std::tolower(chr2);
                                                                   });
        };

        auto case_insensitive_ends_with = [&](std::string const & ext)
        {
            return case_insensitive_string_ends_with(extension, ext);
        };

        return std::ranges::find_if(valid_extensions, case_insensitive_ends_with) != valid_extensions.end();
    }
};

} // namespace raptor

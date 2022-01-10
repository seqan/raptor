// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/io/sequence_file/input.hpp>

#include <raptor/strong_types.hpp>

namespace raptor
{

struct power_of_two_validator
{
    using option_value_type = size_t;

    void operator() (option_value_type const & val) const
    {
        if (!std::has_single_bit(val))
            throw seqan3::validation_error{"The value must be a power of two."};
    }

    static std::string get_help_page_message ()
    {
        return "Value must be a power of two.";
    }
};

class positive_integer_validator
{
public:
    using option_value_type = size_t;

    positive_integer_validator() = default;
    positive_integer_validator(positive_integer_validator const &) = default;
    positive_integer_validator & operator=(positive_integer_validator const &) = default;
    positive_integer_validator(positive_integer_validator &&) = default;
    positive_integer_validator & operator=(positive_integer_validator &&) = default;
    ~positive_integer_validator() = default;

    explicit positive_integer_validator(bool const is_zero_positive_) : is_zero_positive{is_zero_positive_} {}

    void operator()(window const & val) const
    {
        return operator()(val.v);
    }

    void operator()(option_value_type const & val) const
    {
        if (!is_zero_positive && !val)
            throw seqan3::validation_error{"The value must be a positive integer."};
    }

    std::string get_help_page_message () const
    {
        if (is_zero_positive)
            return "Value must be a positive integer or 0.";
        else
            return "Value must be a positive integer.";
    }

private:
    bool is_zero_positive{false};
};

class size_validator
{
public:
    using option_value_type = std::string;

    size_validator() = default;
    size_validator(size_validator const &) = default;
    size_validator & operator=(size_validator const &) = default;
    size_validator(size_validator &&) = default;
    size_validator & operator=(size_validator &&) = default;
    ~size_validator() = default;

    explicit size_validator(std::string const & pattern) : expression{pattern} {}

    void operator()(option_value_type const & cmp) const
    {
        if (!std::regex_match(cmp, expression))
            throw seqan3::validation_error{seqan3::detail::to_string("Value ", cmp, " must be an integer followed by [k,m,g,t] (case insensitive).")};
    }

    template <std::ranges::forward_range range_type>
        requires std::convertible_to<std::ranges::range_value_t<range_type>, option_value_type const &>
    void operator()(range_type const & v) const
    {
         std::for_each(v.begin(), v.end(), [&] (auto cmp) { (*this)(cmp); });
    }

    std::string get_help_page_message() const
    {
        return "Must be an integer followed by [k,m,g,t] (case insensitive).";
    }

private:
    std::regex expression;
};

class bin_validator
{
public:
    using option_value_type = std::vector<std::vector<std::string>>;

    bin_validator() = default;
    bin_validator(bin_validator const &) = default;
    bin_validator & operator=(bin_validator const &) = default;
    bin_validator(bin_validator &&) = default;
    bin_validator & operator=(bin_validator &&) = default;
    ~bin_validator() = default;

    void operator() (option_value_type const & values) const
    {
        if (values.empty())
            throw seqan3::validation_error{"The list of input files cannot be empty."};

        bool const is_minimiser_input = std::filesystem::path{values[0][0]}.extension() == ".minimiser";

        for (std::vector<std::string> const & vector_of_paths : values)
        {
            for (std::string const & value : vector_of_paths)
            {
                std::filesystem::path const file_path{value};

                if (is_minimiser_input && (file_path.extension() != ".minimiser"))
                    throw seqan3::validation_error{"You cannot mix sequence and minimiser files as input."};
                if (std::filesystem::file_size(file_path) == 0u)
                    throw seqan3::validation_error{"The file " + value + " is empty."};

                if (is_minimiser_input)
                    minimiser_file_validator(file_path);
                else
                    sequence_file_validator(file_path);
            }
        }
    }

    std::string get_help_page_message() const
    {
        // Update
        return seqan3::detail::to_string("The file must contain at least one file path per line, with multiple paths "
                                         "being separated by a whitespace. Each line in the file corresponds to one "
                                         "bin. Valid extensions for the paths in the file are [minimiser] when "
                                         " preprocessing, and ", sequence_extensions,
                                         #if defined(SEQAN3_HAS_BZIP2) || defined(SEQAN3_HAS_ZLIB)
                                         ", possibly followed by ", compression_extensions,
                                         #endif
                                         " otherwise. ");
    }

private:
    std::vector<std::string> sequence_extensions{seqan3::detail::valid_file_extensions<typename seqan3::sequence_file_input<>::valid_formats>()};
    std::vector<std::string> compression_extensions{[&] ()
                             {
                                 std::vector<std::string> result;
                                 #ifdef SEQAN3_HAS_BZIP2
                                     result.push_back("bz2");
                                 #endif
                                 #ifdef SEQAN3_HAS_ZLIB
                                     result.push_back("gz");
                                     result.push_back("bgzf");
                                 #endif
                                 return result;
                             }()}; // LCOV_EXCL_LINE
    std::vector<std::string> combined_extensions{[&] ()
                             {
                                 if (compression_extensions.empty())
                                    return sequence_extensions; // LCOV_EXCL_LINE
                                 std::vector<std::string> result;
                                 for (auto && sequence_extension : sequence_extensions)
                                 {
                                     result.push_back(sequence_extension);
                                     for (auto && compression_extension : compression_extensions)
                                         result.push_back(sequence_extension + std::string{'.'} + compression_extension);
                                 }
                                return result;
                             }()};
    seqan3::input_file_validator<> minimiser_file_validator{{"minimiser"}};

public:
    seqan3::input_file_validator<> sequence_file_validator{{combined_extensions}};
};

} // namespace raptor

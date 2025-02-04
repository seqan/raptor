// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides various validator.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <sharg/parser.hpp>

#include <seqan3/io/sequence_file/input.hpp>

#include <raptor/strong_types.hpp>

namespace raptor::detail
{

static inline std::vector<std::string> sequence_extensions{
    seqan3::detail::valid_file_extensions<typename seqan3::sequence_file_input<>::valid_formats>()};

static inline std::vector<std::string> compression_extensions{[]()
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
                                                              }()}; // GCOVR_EXCL_LINE

static inline std::vector<std::string> combined_extensions{
    []()
    {
        std::vector<std::string> result;
        if (compression_extensions.empty())
        {
            result = sequence_extensions; // GCOVR_EXCL_LINE
        }
        else
        {
            for (auto && sequence_extension : sequence_extensions)
            {
                result.push_back(sequence_extension);
                for (auto && compression_extension : compression_extensions)
                {
                    result.push_back(sequence_extension);
                    result.back() += '.';
                    result.back() += compression_extension;
                }
            }
        }
        return result;
    }()};

} // namespace raptor::detail

namespace raptor
{

struct power_of_two_validator
{
    using option_value_type = size_t;

    void operator()(option_value_type const & val) const
    {
        if (!std::has_single_bit(val))
            throw sharg::validation_error{"The value must be a power of two."};
    }

    static std::string get_help_page_message()
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

    explicit positive_integer_validator(bool const is_zero_positive_) : is_zero_positive{is_zero_positive_}
    {}

    void operator()(option_value_type const & val) const
    {
        if (!is_zero_positive && !val)
            throw sharg::validation_error{"The value must be a positive integer."};
    }

    std::string get_help_page_message() const
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

    explicit size_validator(std::string const & pattern) : expression{pattern}
    {}

    void operator()(option_value_type const & cmp) const
    {
        if (!std::regex_match(cmp, expression))
            throw sharg::validation_error{
                seqan3::detail::to_string("Value ",
                                          cmp,
                                          " must be an integer followed by [k,m,g,t] (case insensitive).")};
    }

    template <std::ranges::forward_range range_type>
        requires std::convertible_to<std::ranges::range_value_t<range_type>, option_value_type const &>
    void operator()(range_type const & v) const
    {
        std::for_each(v.begin(),
                      v.end(),
                      [&](auto cmp)
                      {
                          (*this)(cmp);
                      });
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

    void operator()(option_value_type const & values) const
    {
        if (values.empty())
            throw sharg::validation_error{"The list of input files cannot be empty."};

        bool const is_minimiser_input = std::filesystem::path{values[0][0]}.extension() == ".minimiser";

        for (std::vector<std::string> const & vector_of_paths : values)
        {
            for (std::string const & value : vector_of_paths)
            {
                std::filesystem::path const file_path{value};

                if (is_minimiser_input && (file_path.extension() != ".minimiser"))
                    throw sharg::validation_error{"You cannot mix sequence and minimiser files as input."};
                if (!std::filesystem::exists(file_path))
                    throw sharg::validation_error{"The file " + value + " does not exist."};
                if (std::filesystem::file_size(file_path) == 0u)
                    throw sharg::validation_error{"The file " + value + " is empty."};

                if (is_minimiser_input)
                    minimiser_file_validator(file_path);
                else
                    sequence_file_validator(file_path);
            }
        }
    }

    std::string get_help_page_message() const
    {
        return seqan3::detail::to_string("The file must contain at least one file path per line, with multiple paths "
                                         "being separated by a whitespace. Each line in the file corresponds to one "
                                         "bin. Valid extensions for the paths in the file are [minimiser] when "
                                         " using preprocessed input from \\fBraptor prepare\\fP, and ",
                                         raptor::detail::sequence_extensions,
#if defined(SEQAN3_HAS_BZIP2) || defined(SEQAN3_HAS_ZLIB)
                                         ", possibly followed by ",
                                         raptor::detail::compression_extensions,
#endif
                                         ". ");
    }

private:
    sharg::input_file_validator minimiser_file_validator{{"minimiser"}};

public:
    sharg::input_file_validator sequence_file_validator{raptor::detail::combined_extensions};
};

class output_directory_validator
{
public:
    using option_value_type = std::string;

    output_directory_validator() = default;
    output_directory_validator(output_directory_validator const &) = default;
    output_directory_validator & operator=(output_directory_validator const &) = default;
    output_directory_validator(output_directory_validator &&) = default;
    output_directory_validator & operator=(output_directory_validator &&) = default;
    ~output_directory_validator() = default;

    void operator()(option_value_type const & value) const
    {
        std::filesystem::path const out_dir{value};
        std::error_code ec{};
        std::filesystem::create_directories(out_dir, ec);
        if (ec)
            // GCOVR_EXCL_START
            throw sharg::validation_error{
                sharg::detail::to_string("Failed to create directory\"", out_dir.c_str(), "\": ", ec.message())};
        // GCOVR_EXCL_STOP

        validator(out_dir);
    }

    std::string get_help_page_message() const
    {
        return "A valid path for the output directory.";
    }

private:
    sharg::output_directory_validator validator{};
};

class output_file_validator
{
public:
    using option_value_type = std::string;

    output_file_validator() = default;
    output_file_validator(output_file_validator const &) = default;
    output_file_validator & operator=(output_file_validator const &) = default;
    output_file_validator(output_file_validator &&) = default;
    output_file_validator & operator=(output_file_validator &&) = default;
    ~output_file_validator() = default;

    void operator()(option_value_type const & value) const
    {
        std::filesystem::path const out_path{value};
        std::filesystem::path const out_dir{out_path.parent_path()};
        if (!out_dir.empty())
        {
            // GCOVR_EXCL_START
            std::error_code ec{};
            std::filesystem::create_directories(out_dir, ec);
            if (ec)
                throw sharg::validation_error{
                    sharg::detail::to_string("Failed to create directory \"", out_dir.c_str(), "\": ", ec.message())};
            // GCOVR_EXCL_STOP
        }

        validator(out_path);
    }

    std::string get_help_page_message() const
    {
        return "A valid path for the output file. Write permissions must be granted.";
    }

private:
    sharg::output_file_validator validator{sharg::output_file_open_options::open_or_create};
};

class sequence_file_validator : public sharg::input_file_validator
{
private:
    using base_t = sharg::input_file_validator;

public:
    using base_t::base_t;

    std::string get_help_page_message() const
    {
        return seqan3::detail::to_string(
            "The input file must exist and read permissions must be granted. Valid file extensions are ",
            raptor::detail::sequence_extensions,
#if defined(SEQAN3_HAS_BZIP2) || defined(SEQAN3_HAS_ZLIB)
            ", possibly followed by ",
            raptor::detail::compression_extensions,
#endif
            ". ");
    }
};

} // namespace raptor

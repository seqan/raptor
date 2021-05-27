#pragma once

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/input.hpp>

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
    using option_value_type = std::vector<std::filesystem::path>;

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

        for (auto && value : values)
        {
            try
            {
                sequence_file_validator(value);
            }
            catch (seqan3::validation_error const & exception)
            {
                if (value.extension() == ".minimiser")
                    minimiser_file_validator(value);
                else if (values.size() == 1u)
                {
                    std::ifstream list_of_files{value};
                    std::string line;
                    while (std::getline(list_of_files, line))
                    {
                        if (!line.empty())
                        {
                            std::filesystem::path bin_path{line};
                            if (bin_path.extension() == ".minimiser")
                                minimiser_file_validator(bin_path);
                            else
                                sequence_file_validator(bin_path);
                        }
                    }
                }
                else
                    throw exception;
            }
        }

        bool const is_minimiser_input = values[0].extension() == ".minimiser";

        for (auto && value : values)
            if (is_minimiser_input != (value.extension() == ".minimiser"))
                throw seqan3::validation_error{"You cannot mix sequence and minimiser files as input."};
    }

    std::string get_help_page_message() const
    {
        return seqan3::detail::to_string("The input file must exist and read permissions must be granted. Valid file "
                                         "extensions for bin files are: [minimiser], or ", sequence_extensions,
                                         #if defined(SEQAN3_HAS_BZIP2) || defined(SEQAN3_HAS_ZLIB)
                                         " possibly followed by: ", compression_extensions, ". ",
                                         #else
                                         ". ",
                                         #endif
                                         "All other extensions will be assumed to contain one line per path to a bin.");
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
                             }()};
    std::vector<std::string> combined_extensions{[&] ()
                             {
                                 if (compression_extensions.empty())
                                    return sequence_extensions;
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

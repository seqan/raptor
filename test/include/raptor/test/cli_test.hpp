// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <gtest/gtest.h>

#include <charconv>
#include <fstream>

#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/utility/views/repeat_n.hpp>

#include <hibf/contrib/std/zip_view.hpp>

#include <raptor/index.hpp>

#ifndef RAPTOR_ASSERT_ZERO_EXIT
#    define RAPTOR_ASSERT_ZERO_EXIT(arg) ASSERT_EQ(arg.exit_code, 0) << "Command: " << arg.command
#endif
#ifndef RAPTOR_ASSERT_FAIL_EXIT
#    define RAPTOR_ASSERT_FAIL_EXIT(arg) ASSERT_NE(arg.exit_code, 0) << "Command: " << arg.command
#endif

// Provides functions for CLI test implementation.
struct cli_test : public ::testing::Test
{
private:
    // Holds the original work directory where Gtest has been started.
    std::filesystem::path original_workdir{};

protected:
    // Result struct for captured streams and exit code.
    struct cli_test_result
    {
        std::string out{};
        std::string err{};
        std::string command{};
        int exit_code{};
    };

    // Invoke the app execution. The command line call should be given as separate parameters.
    template <typename... CommandItemTypes>
    cli_test_result execute_app(CommandItemTypes &&... command_items)
    {
        cli_test_result result{};

        // Assemble the command string and disable version check.
        result.command = [&command_items...]()
        {
            std::ostringstream command{};
            command << "SHARG_NO_VERSION_CHECK=1 " << BINDIR;
            (void)((command << command_items << ' '), ...); // (void) silences "state has no effect" warning if no items
            return command.str();
        }();

        // Always capture the output streams.
        testing::internal::CaptureStdout();
        testing::internal::CaptureStderr();

        // Run the command and return results.
        result.exit_code = std::system(result.command.c_str());
        result.out = testing::internal::GetCapturedStdout();
        result.err = testing::internal::GetCapturedStderr();
        return result;
    }

    // Generate the full path of a test input file that is provided in the data directory.
    static std::filesystem::path data(std::string const & filename)
    {
        return std::filesystem::path{std::string{DATADIR}}.concat(filename);
    }

    // Create an individual work directory for the current test.
    void SetUp() override
    {
        // Assemble the directory name.
        ::testing::TestInfo const * const info = ::testing::UnitTest::GetInstance()->current_test_info();
        std::filesystem::path const test_dir{std::string{OUTPUTDIR} + std::string{info->test_case_name()}
                                             + std::string{"."} + std::string{info->name()}};
        try
        {
            std::filesystem::remove_all(test_dir);              // delete the directory if it exists
            std::filesystem::create_directories(test_dir);      // create the new empty directory
            original_workdir = std::filesystem::current_path(); // store original work dir path
            std::filesystem::current_path(test_dir);            // change the work dir
        }
        catch (std::exception const & exc)
        {
            FAIL() << "Failed to set up the test directory " << test_dir << ":\n" << exc.what();
        }
    }

    // Switch back to the initial work directory.
    void TearDown() override
    {
        try
        {
            std::filesystem::current_path(original_workdir); // restore the original work dir
        }
        catch (std::exception const & exc)
        {
            FAIL() << "Failed to set the work directory to " << original_workdir << ":\n" << exc.what();
        }
    }
};

struct raptor_base : public cli_test
{
    struct strong_bool
    {
        enum
        {
            yes,
            no
        } value;

        strong_bool(decltype(value) v) : value(v)
        {}

        explicit operator bool() const
        {
            return value == strong_bool::yes;
        }

        bool operator!() const
        {
            return value == strong_bool::no;
        }
    };

    struct is_hibf : strong_bool
    {
        using strong_bool::strong_bool;
        using strong_bool::value;
    };

    struct compare_extension : strong_bool
    {
        using strong_bool::strong_bool;
        using strong_bool::value;
    };

    struct is_empty : strong_bool
    {
        using strong_bool::strong_bool;
        using strong_bool::value;
    };

    struct is_preprocessed : strong_bool
    {
        using strong_bool::strong_bool;
        using strong_bool::value;
    };

    static inline auto get_repeated_bins(size_t const repetitions) noexcept
    {
        using vec_t = std::vector<std::string>;

        if (repetitions == 0)
            return seqan3::views::repeat_n(vec_t{cli_test::data("bin1.fa").string()}, 1) | std::views::join;

        return seqan3::views::repeat_n(vec_t{cli_test::data("bin1.fa"),
                                             cli_test::data("bin2.fa"),
                                             cli_test::data("bin3.fa"),
                                             cli_test::data("bin4.fa")},
                                       repetitions)
             | std::views::join;
    }

    static inline std::filesystem::path const
    ibf_path(size_t const number_of_repetitions, size_t const window_size, is_hibf const hibf = is_hibf::no) noexcept
    {
        std::string name{};
        name += std::to_string(std::max<int>(1, number_of_repetitions * 4));
        name += "bins";
        name += std::to_string(window_size);
        name += "window.";
        name += hibf ? "hibf" : "index";
        return cli_test::data(name);
    }

    static inline std::filesystem::path const pack_path(size_t const number_of_repetitions) noexcept
    {
        std::string name{};
        name += std::to_string(std::max<int>(1, number_of_repetitions * 4));
        name += "bins.pack";
        return cli_test::data(name);
    }

    static inline std::string const string_from_file(std::filesystem::path const & path,
                                                     std::ios_base::openmode const mode = std::ios_base::in)
    {
        std::ifstream file_stream(path, mode);
        if (!file_stream.is_open())
            throw std::logic_error{"Cannot open " + path.string()};
        std::stringstream file_buffer;
        file_buffer << file_stream.rdbuf();
        return {file_buffer.str()};
    }

    // Good example for printing tables: https://en.cppreference.com/w/cpp/io/ios_base/width
    static inline std::string const debug_ibfs(seqan::hibf::interleaved_bloom_filter const & expected_ibf,
                                               seqan::hibf::interleaved_bloom_filter const & actual_ibf)
    {
        std::stringstream result{};
        result << ">>>IBFs differ<<<\n";
        result.setf(std::ios::left, std::ios::adjustfield);

        result.width(22);
        result << "#Member accessor";
        result.width(15);
        result << "Expected value";
        result.width(13);
        result << "Actual value";
        result << '\n';

        result.width(22);
        result << "bin_count()";
        result.width(15);
        result << expected_ibf.bin_count();
        result.width(13);
        result << actual_ibf.bin_count();
        result << '\n';

        result.width(22);
        result << "bin_size()";
        result.width(15);
        result << expected_ibf.bin_size();
        result.width(13);
        result << actual_ibf.bin_size();
        result << '\n';

        result.width(22);
        result << "hash_function_count()";
        result.width(15);
        result << expected_ibf.hash_function_count();
        result.width(13);
        result << actual_ibf.hash_function_count();
        result << '\n';

        result.width(22);
        result << "bit_size()";
        result.width(15);
        result << expected_ibf.bit_size();
        result.width(13);
        result << actual_ibf.bit_size();
        result << '\n';

        return result.str();
    }

    template <typename data_t = raptor::index_structure::ibf>
    static inline void compare_index(std::filesystem::path const & expected_result,
                                     std::filesystem::path const & actual_result,
                                     compare_extension const compare_ext = compare_extension::yes)
    {
        constexpr bool is_ibf = std::same_as<data_t, raptor::index_structure::ibf>;
        constexpr bool is_hibf = std::same_as<data_t, raptor::index_structure::hibf>;

        static_assert(is_ibf || is_hibf);

        raptor::raptor_index<data_t> expected_index{}, actual_index{};

        {
            std::ifstream is{expected_result, std::ios::binary};
            cereal::BinaryInputArchive iarchive{is};
            iarchive(expected_index);
        }
        {
            std::ifstream is{actual_result, std::ios::binary};
            cereal::BinaryInputArchive iarchive{is};
            iarchive(actual_index);
        }

        EXPECT_EQ(expected_index.window_size(), actual_index.window_size());
        EXPECT_EQ(expected_index.shape(), actual_index.shape());
        EXPECT_EQ(expected_index.parts(), actual_index.parts());

        if constexpr (is_ibf)
        {
            auto const &expected_ibf{expected_index.ibf()}, actual_ibf{actual_index.ibf()};
            EXPECT_TRUE(expected_ibf == actual_ibf) << debug_ibfs(expected_ibf, actual_ibf);
        }
        else
        {
            auto const &expected_ibfs{expected_index.ibf().ibf_vector}, actual_ibfs{actual_index.ibf().ibf_vector};
            for (auto const & expected_ibf : expected_ibfs)
            {
                ASSERT_TRUE(std::ranges::find(actual_ibfs, expected_ibf) != actual_ibfs.end());
            }
        }

        auto const &all_expected_bins{expected_index.bin_path()}, all_actual_bins{actual_index.bin_path()};
        EXPECT_EQ(std::ranges::distance(all_expected_bins), std::ranges::distance(all_actual_bins));

        if constexpr (is_ibf)
        {
            for (auto const && [expected_list, actual_list] :
                 seqan::stl::views::zip(all_expected_bins, all_actual_bins))
            {
                EXPECT_TRUE(std::ranges::distance(expected_list) > 0);
                for (auto const && [expected_file, actual_file] : seqan::stl::views::zip(expected_list, actual_list))
                {
                    std::filesystem::path const expected_path(expected_file);
                    std::filesystem::path const actual_path(actual_file);
                    if (compare_ext)
                        ASSERT_EQ(expected_path.filename(), actual_path.filename());
                    else
                        ASSERT_EQ(expected_path.stem(), actual_path.stem());
                }
            }
        }
        else
        {
            auto filenames = std::views::transform(
                [compare_ext](std::vector<std::string> const & filename_list)
                {
                    std::vector<std::string> result{};
                    if (compare_ext)
                        for (auto const & filename : filename_list)
                            result.emplace_back(std::filesystem::path{filename}.filename().string());
                    else
                        for (auto const & filename : filename_list)
                            result.emplace_back(std::filesystem::path{filename}.stem().string());
                    return result;
                });

            auto expected_filenames_view = all_expected_bins | filenames;
            auto actual_filenames_view = all_actual_bins | filenames;

            std::vector<std::vector<std::string>> expected_filenames(expected_filenames_view.begin(),
                                                                     expected_filenames_view.end());
            std::vector<std::vector<std::string>> actual_filenames(actual_filenames_view.begin(),
                                                                   actual_filenames_view.end());
            std::ranges::sort(expected_filenames);
            std::ranges::sort(actual_filenames);
            EXPECT_RANGE_EQ(expected_filenames, actual_filenames);
        }
    }

    static inline void compare_search(size_t const number_of_repeated_bins,
                                      size_t const number_of_errors,
                                      std::string_view const filename,
                                      is_empty const empty = is_empty::no,
                                      is_preprocessed const preprocessed = is_preprocessed::no)
    {
        size_t const number_of_bins{std::max<size_t>(1u, number_of_repeated_bins * 4u)};
        std::string_view const missed_bin = number_of_errors ? "none" : preprocessed ? "bin4.minimiser" : "bin4.fa";
        std::ifstream search_result{filename.data()};
        std::string line;
        std::vector<uint64_t> expected_hits;
        std::vector<uint64_t> actual_hits;
        uint64_t tmp{};

        // Skip parameter information
        while (std::getline(search_result, line) && line.starts_with("##"))
        {}

        // Parse user bin IDs, skipping missed_bin
        for (size_t i = 0; i < number_of_bins; ++i)
        {
            std::string_view line_view{line};
            if (!empty && !line_view.ends_with(missed_bin))
            {
                std::from_chars(line_view.data() + 1u, line_view.data() + line_view.find('\t'), tmp);
                expected_hits.push_back(tmp);
            }
            ASSERT_TRUE(std::getline(search_result, line));
        }

        ASSERT_EQ(line, "#QUERY_NAME\tUSER_BINS");

        std::string const query_prefix{"query"};
        for (char i : {'3', '1', '2'})
        {
            ASSERT_TRUE(std::getline(search_result, line));
            std::string_view line_view{line};

            ASSERT_EQ(line_view.substr(0, query_prefix.size() + 2u), query_prefix + i + '\t');

            actual_hits.clear();

            for (auto && hit : std::views::split(line_view.substr(query_prefix.size() + 2u), ','))
            {
// GCC 11 does not implement P2210R2, hence GCC 11's split_view is actually a lazy_split_view.
#if defined(__GNUC__) && !defined(__llvm__) && !defined(__INTEL_COMPILER) && (__GNUC__ > 10) && (__GNUC__ < 12)
                std::stringstream buf;
                for (auto const & inner_view : hit)
                {
                    buf << inner_view;
                }
                auto view = buf.view();
                std::from_chars(view.data(), view.data() + view.size(), tmp);
#else
                std::from_chars(hit.data(), hit.data() + hit.size(), tmp);
#endif
                actual_hits.push_back(tmp);
            }
            std::ranges::sort(actual_hits);
            ASSERT_EQ(expected_hits, actual_hits);
        }

        EXPECT_FALSE(std::getline(search_result, line));
    }
};

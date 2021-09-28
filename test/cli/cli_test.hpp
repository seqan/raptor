// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <gtest/gtest.h>

#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/utility/views/zip.hpp>

#include <raptor/index.hpp>

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
        int exit_code{};
    };

    // Invoke the app execution. The command line call should be given as separate parameters.
    template <typename... CommandItemTypes>
    cli_test_result execute_app(CommandItemTypes &&... command_items)
    {
        cli_test_result result{};

        // Assemble the command string and disable version check.
        std::ostringstream command{};
        command << "SEQAN3_NO_VERSION_CHECK=1 " << BINDIR;
        int a[] = {0, ((void)(command << command_items << ' '), 0) ... };
        (void) a;

        // Always capture the output streams.
        testing::internal::CaptureStdout();
        testing::internal::CaptureStderr();

        // Run the command and return results.
        result.exit_code = std::system(command.str().c_str());
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
        std::filesystem::path const test_dir{std::string{OUTPUTDIR} +
                                             std::string{info->test_case_name()} +
                                             std::string{"."} +
                                             std::string{info->name()}};
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
            std::filesystem::current_path(original_workdir);    // restore the original work dir
        }
        catch (std::exception const & exc)
        {
            FAIL() << "Failed to set the work directory to " << original_workdir << ":\n" << exc.what();
        }
    }
};

struct raptor_base : public cli_test
{
    static inline std::string const repeat_bins(size_t const repetitions) noexcept
    {
        if (repetitions == 0)
            return cli_test::data("bin1.fa").string();

        std::string result{};

        for (size_t i{0}; i < repetitions; ++i)
        {
            result += cli_test::data("bin1.fa");
            result += ' ';
            result += cli_test::data("bin2.fa");
            result += ' ';
            result += cli_test::data("bin3.fa");
            result += ' ';
            result += cli_test::data("bin4.fa");
            result += ' ';
        }

        return result;
    }

    static inline std::filesystem::path const ibf_path(size_t const number_of_repetitions,
                                                       size_t const window_size,
                                                       bool const compressed = false,
                                                       bool const hibf = false) noexcept
    {
        std::string name{};
        name += std::to_string(std::max<int>(1, number_of_repetitions * 4));
        name += "bins";
        name += std::to_string(window_size);
        name += compressed ? "windowc." : "window.";
        name += hibf ? "hibf" : "index";
        return cli_test::data(name);
    }

    static inline std::filesystem::path const search_result_path(size_t const number_of_repetitions,
                                                                 size_t const window_size,
                                                                 size_t const number_of_errors,
                                                                 bool const socks = false,
                                                                 bool const empty = false,
                                                                 bool const hibf = false) noexcept
    {
        std::string name{};
        name += std::to_string(std::max<int>(1, number_of_repetitions * 4));
        name += "bins";
        if (empty)
        {
            name += "empty";
        }
        else
        {
            name += std::to_string(window_size);
            name += "window";
            name += std::to_string(number_of_errors);
            name += "error";
        }
        if (socks)
            name += "socks.out";
        else if (hibf)
            name += "hibf.out";
        else
            name += ".out";
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
    template <seqan3::data_layout layout = seqan3::data_layout::uncompressed>
    static inline std::string const debug_ibfs(seqan3::interleaved_bloom_filter<layout> const & expected_ibf,
                                               seqan3::interleaved_bloom_filter<layout> const & actual_ibf)
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
    static inline void compare_results(std::filesystem::path const & expected_result,
                                       std::filesystem::path const & actual_result,
                                       bool const compare_extension = true)
    {
        constexpr bool is_ibf = std::same_as<data_t, raptor::index_structure::ibf> ||
                                std::same_as<data_t, raptor::index_structure::ibf_compressed>;
        constexpr bool is_hibf = std::same_as<data_t, raptor::index_structure::hibf> ||
                                 std::same_as<data_t, raptor::index_structure::hibf_compressed>;

        static_assert(is_ibf || is_hibf);

        constexpr auto data_layout = [is_ibf] () constexpr
        {
            if constexpr (is_ibf)
                return data_t::data_layout_mode;
            else
                return data_t::ibf_t::data_layout_mode;
        }();

        raptor::raptor_index<data_t> expected_index, actual_index{};

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
        EXPECT_EQ(expected_index.compressed(), actual_index.compressed());

        if constexpr(is_ibf)
        {
            auto const & expected_ibf{expected_index.ibf()}, actual_ibf{actual_index.ibf()};
            EXPECT_TRUE(expected_ibf == actual_ibf) << debug_ibfs<data_layout>(expected_ibf, actual_ibf);
        }
        else
        {
            auto const & expected_hibf{expected_index.ibf()}, actual_hibf{actual_index.ibf()};
            for (auto const && [expected_ibf, actual_ibf] : seqan3::views::zip(expected_hibf, actual_hibf))
            {
                EXPECT_TRUE(expected_ibf == actual_ibf) << debug_ibfs<data_layout>(expected_ibf, actual_ibf);
            }
        }

        auto const & all_expected_bins{expected_index.bin_path()}, all_actual_bins{actual_index.bin_path()};
        EXPECT_EQ(std::ranges::distance(all_expected_bins), std::ranges::distance(all_actual_bins));

        for (auto const && [expected_list, actual_list] : seqan3::views::zip(all_expected_bins, all_actual_bins))
        {
            EXPECT_TRUE(std::ranges::distance(expected_list) > 0);
            for (auto const && [expected_file, actual_file] : seqan3::views::zip(expected_list, actual_list))
            {
                std::filesystem::path const expected_path(expected_file);
                std::filesystem::path const actual_path(actual_file);
                if (compare_extension)
                    EXPECT_EQ(expected_path.filename(), actual_path.filename());
                else
                    EXPECT_EQ(expected_path.stem(), actual_path.stem());
            }
        }
    }
};

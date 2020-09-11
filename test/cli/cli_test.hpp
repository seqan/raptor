#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <cstdlib>               // system calls
#include <seqan3/std/filesystem> // test directory creation
#include <sstream>               // ostringstream
#include <string>                // strings

// Include the EXPECT_RANGE_EQ macro for better information if range elements differ.
#include <seqan3/test/expect_range_eq.hpp>

#pragma once

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

struct raptor : public cli_test
{
    std::string expand_bins(size_t const bins)
    {
        std::string result{};

        for (size_t i{0}; i < bins; ++i)
        {
            std::string filename = "example_data/" +
                                   std::to_string(bins) +
                                   "/bins" +
                                   "/bin_" +
                                   std::string(std::to_string(bins).size() - std::to_string(i).size(), '0') + // prepend zeros
                                   std::to_string(i) +
                                   ".fasta";
            result += cli_test::data(filename);
            result += ' ';
        }

        return result;
    }
};

struct raptor_build : public raptor, public testing::WithParamInterface<std::tuple<size_t, bool>> {};
struct raptor_search : public raptor, public testing::WithParamInterface<std::tuple<size_t, size_t, size_t>> {};

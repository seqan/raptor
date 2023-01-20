// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/test/file_access.hpp>

#include <raptor/build/tmp_directory.hpp>

TEST(tmp_directory, unique)
{
    raptor::tmp_directory t1{};
    raptor::tmp_directory t2{};

    EXPECT_NE(t1.path(), t2.path());

    EXPECT_TRUE(std::filesystem::exists(t1.path()));
    EXPECT_TRUE(std::filesystem::exists(t2.path()));

    // check if created folders are empty
    EXPECT_TRUE(t1.empty());
    EXPECT_TRUE(t2.empty());

    // checking they are inside of /tmp
    EXPECT_TRUE(std::filesystem::equivalent(std::filesystem::temp_directory_path(), t1.path().parent_path()));
    EXPECT_TRUE(std::filesystem::equivalent(std::filesystem::temp_directory_path(), t2.path().parent_path()));
}

TEST(tmp_directory, move_constructible)
{
    raptor::tmp_directory t1{};
    raptor::tmp_directory t2{};
    raptor::tmp_directory t3{std::move(t2)};

    EXPECT_TRUE(std::filesystem::exists(t1.path()));
    EXPECT_TRUE(std::filesystem::exists(t3.path()));
    EXPECT_TRUE(t1.empty());
    EXPECT_TRUE(t3.empty());

    EXPECT_NE(t1.path(), t3.path());
    raptor::tmp_directory t4(std::move(t1));

    EXPECT_TRUE(std::filesystem::exists(t4.path()));
    EXPECT_NE(t3.path(), t4.path());
}

TEST(tmp_directory, move_assignable)
{
    std::filesystem::path p1{};
    std::filesystem::path p2{};
    std::filesystem::path p3{};

    {
        raptor::tmp_directory t1{};
        raptor::tmp_directory t2{};
        raptor::tmp_directory t3{};

        p1 = t1.path();
        p2 = t2.path();
        p3 = t3.path();

        t3 = std::move(t2);

        EXPECT_NE(t1.path(), t3.path());

        EXPECT_TRUE(std::filesystem::exists(t1.path()));
        EXPECT_TRUE(std::filesystem::exists(t3.path()));
    }

    EXPECT_FALSE(std::filesystem::exists(p1));
    EXPECT_FALSE(std::filesystem::exists(p2));
    EXPECT_FALSE(std::filesystem::exists(p3));
}

// check destructor does all its cleanups
TEST(tmp_directory, cleanup_on_destruction)
{
    std::filesystem::path path{};
    {
        raptor::tmp_directory t1{};
        path = t1.path();

        {
            std::ofstream os{path / "file1", std::ios::out};
            os << "some data";
        }

        {
            std::filesystem::create_directory(path / "somefolder");
            std::ofstream os{path / "somefolder/file2"};
            os << "other data";
        }

        EXPECT_FALSE(t1.empty());

        EXPECT_TRUE(std::filesystem::exists(path));
        EXPECT_TRUE(std::filesystem::exists(path / "file1"));
        EXPECT_TRUE(std::filesystem::exists(path / "somefolder/file2"));
    }

    EXPECT_FALSE(std::filesystem::exists(path));
    EXPECT_FALSE(std::filesystem::exists(path / "file1"));
    EXPECT_FALSE(std::filesystem::exists(path / "somefolder/file2"));
}

TEST(tmp_directory_throw, directory_not_writeable)
{
    // create a temporary folder that will mimic the normal tmp folder
    raptor::tmp_directory temporary_tmp_folder;
    setenv("TMPDIR", temporary_tmp_folder.path().c_str(), 1); // name, value, overwrite

    // make temporary_tmp_folder read only
    std::filesystem::permissions(temporary_tmp_folder.path(),
                                 std::filesystem::perms::owner_write,
                                 std::filesystem::perm_options::remove);

    // The actual test
    if (!seqan3::test::write_access(temporary_tmp_folder.path())) // Do not execute with root permissions.
    {
        EXPECT_THROW(raptor::tmp_directory{}, std::filesystem::filesystem_error);
    }

    // give temporary_tmp_folder write permissions back
    std::filesystem::permissions(temporary_tmp_folder.path(),
                                 std::filesystem::perms::owner_write,
                                 std::filesystem::perm_options::add);
}

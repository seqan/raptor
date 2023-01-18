// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#if defined(__APPLE__)
#    include <unistd.h>
#elif defined(_WIN32)
#    include <cstring>
#    include <io.h>
#else
#    include <cstdlib>
#endif

#include <cassert>
#include <filesystem>
#include <optional>

#if defined(_WIN32)
namespace
{

char * mkdtemp(char * template_name)
{
    if (_mktemp_s(template_name, strlen(template_name) + 1))
        return nullptr;

    if (std::filesystem::create_directories(template_name))
        return template_name;

    return nullptr;
}

} // namespace
#endif

namespace raptor
{

class tmp_directory
{
public:
    tmp_directory()
    {
        std::string tmp_template{(std::filesystem::temp_directory_path() / "raptor_XXXXXXXX").string()};
        char * const path_str = tmp_template.data();

        if (char * f = mkdtemp(path_str); f == nullptr)
        {
            throw std::filesystem::filesystem_error("Could not create temporary directory with mkdtemp!",
                                                    path_str,
                                                    std::make_error_code(std::errc::bad_file_descriptor));
        }
        directory_path = path_str;
    }

    tmp_directory(tmp_directory const &) = delete;
    tmp_directory & operator=(tmp_directory const &) = delete;

    tmp_directory(tmp_directory && other) : directory_path{std::exchange(other.directory_path, std::nullopt)}
    {}

    tmp_directory & operator=(tmp_directory && other)
    {
        clean();
        directory_path = std::exchange(other.directory_path, std::nullopt);
        return *this;
    }

    ~tmp_directory() noexcept(false)
    {
        clean();
    }

    std::filesystem::path path() const
    {
        assert(directory_path);

        return directory_path.value();
    }

    bool empty() const
    {
        assert(directory_path);

        return std::filesystem::exists(directory_path.value()) && std::filesystem::is_empty(directory_path.value());
    }

private:
    void clean()
    {
        if (!directory_path)
            return;

        std::filesystem::remove_all(directory_path.value());
        directory_path = std::nullopt;
    }

    std::optional<std::filesystem::path> directory_path;
};

} // namespace raptor

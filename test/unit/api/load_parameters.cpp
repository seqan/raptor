// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <string>

#include <raptor/index.hpp>

using hibf_index_t = raptor::raptor_index<raptor::index_structure::hibf>;

static std::filesystem::path data(std::string const & filename)
{
    return std::filesystem::path{std::string{DATADIR}}.concat(filename);
}

template <typename index_t>
static index_t load(std::filesystem::path const & path, bool const parameters_only)
{
    index_t index{};
    std::ifstream stream{path, std::ios::binary};
    cereal::BinaryInputArchive archive{stream};
    if (parameters_only)
        index.load_parameters(archive);
    else
        archive(index);
    return index;
}

// load_parameters must yield the same parameters as loading the whole index.
static void expect_same_parameters(hibf_index_t const & parameters, hibf_index_t const & full)
{
    EXPECT_EQ(parameters.window_size(), full.window_size());
    EXPECT_EQ(parameters.shape(), full.shape());
    EXPECT_EQ(parameters.parts(), full.parts());
    EXPECT_EQ(parameters.bin_path(), full.bin_path());
    EXPECT_EQ(parameters.fpr(), full.fpr());
    EXPECT_EQ(parameters.is_hibf(), full.is_hibf());
    EXPECT_TRUE(parameters.config() == full.config());
    EXPECT_EQ(parameters.original_number_of_ibfs(), full.original_number_of_ibfs());
    EXPECT_EQ(parameters.number_of_resized_ibfs(), full.number_of_resized_ibfs());
    ASSERT_EQ(parameters.number_of_ibfs(), full.number_of_ibfs());
    for (size_t ibf_idx = 0; ibf_idx < full.number_of_ibfs(); ++ibf_idx)
        EXPECT_EQ(parameters.was_resized(ibf_idx), full.was_resized(ibf_idx)) << "IBF " << ibf_idx;
}

TEST(load_parameters, hibf)
{
    std::filesystem::path const path = data("three_levels.hibf");
    hibf_index_t const parameters = load<hibf_index_t>(path, true);
    hibf_index_t const full = load<hibf_index_t>(path, false);

    expect_same_parameters(parameters, full);

    // The data itself is not loaded.
    EXPECT_TRUE(parameters.ibf().ibf_vector.empty());
    // The resize bookkeeping is as large as the HIBF.
    EXPECT_EQ(full.number_of_ibfs(), full.ibf().ibf_vector.size());
    EXPECT_EQ(parameters.number_of_ibfs(), 49u);
    EXPECT_EQ(parameters.original_number_of_ibfs(), 49u);
    EXPECT_EQ(parameters.number_of_resized_ibfs(), 0u);
}

// The resize bookkeeping is serialised for an IBF index as well, but empty.
TEST(load_parameters, ibf)
{
    std::filesystem::path const path = data("64bins19window.index");
    hibf_index_t const parameters = load<hibf_index_t>(path, true);
    auto const full = load<raptor::raptor_index<raptor::index_structure::ibf>>(path, false);

    EXPECT_EQ(parameters.window_size(), full.window_size());
    EXPECT_EQ(parameters.shape(), full.shape());
    EXPECT_EQ(parameters.parts(), full.parts());
    EXPECT_EQ(parameters.bin_path(), full.bin_path());
    EXPECT_EQ(parameters.fpr(), full.fpr());
    EXPECT_FALSE(parameters.is_hibf());
    EXPECT_EQ(parameters.number_of_ibfs(), 0u);
    EXPECT_EQ(parameters.original_number_of_ibfs(), 0u);
    EXPECT_EQ(parameters.number_of_resized_ibfs(), 0u);
}

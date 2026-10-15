// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "inspect_test.hpp"

struct inspect_ibf : public inspect_test, public testing::WithParamInterface<std::tuple<size_t, size_t>>
{};

// The shipped IBF indices are built with `raptor build --kmer 19 --window <W>` and default parameters otherwise.
TEST_P(inspect_ibf, text)
{
    auto const [number_of_repeated_bins, window_size] = GetParam();
    std::filesystem::path const index = ibf_path(number_of_repeated_bins, window_size);

    cli_test_result const result = execute_app("raptor", "inspect", "--index", index);
    std::string const expected = "### Index\n"
                                 "## Index = \""
                               + index.string()
                               + "\"\n"
                                 "## Index version = 3\n"
                                 "## Index type = IBF\n"
                                 "## Index parts = 1\n"
                                 "## Index size [bytes] = "
                               + std::to_string(std::filesystem::file_size(index))
                               + "\n"
                                 "### Minimiser parameters\n"
                                 "## Window size = "
                               + std::to_string(window_size)
                               + "\n"
                                 "## Shape = 1111111111111111111\n"
                                 "## Shape size (length) = 19\n"
                                 "## Shape count (number of 1s) = 19\n"
                                 "### Bloom filter parameters\n"
                                 "## False positive rate = 0.05\n"
                                 "### User bins\n"
                                 "## Number of user bins = "
                               + std::to_string(std::max<size_t>(1u, number_of_repeated_bins * 4u)) + "\n";
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);
}

TEST_P(inspect_ibf, user_bins)
{
    auto const [number_of_repeated_bins, window_size] = GetParam();
    std::filesystem::path const index = ibf_path(number_of_repeated_bins, window_size);

    cli_test_result const result = execute_app("raptor", "inspect", "--index", index, "--user-bins");
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);

    text_output const output = parse_text(result.out);
    EXPECT_EQ(output.user_bins, user_bin_lines(load_index<ibf_index_t>(index).bin_path()));
    EXPECT_TRUE(output.ibfs.empty());
}

TEST_P(inspect_ibf, detailed)
{
    auto const [number_of_repeated_bins, window_size] = GetParam();
    std::filesystem::path const index = ibf_path(number_of_repeated_bins, window_size);

    cli_test_result const result = execute_app("raptor", "inspect", "--index", index, "--detailed");
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);

    ibf_index_t const loaded = load_index<ibf_index_t>(index);
    auto const & ibf = loaded.ibf();
    expected_ibf const expected = expected_from(ibf, 0u);

    text_output const output = parse_text(result.out);
    EXPECT_EQ(output.sections.back(), "Structure");
    EXPECT_EQ(output.values.at("Levels"), "1");
    EXPECT_EQ(output.values.at("Technical bins"), std::to_string(ibf.bin_count()));
    EXPECT_EQ(output.values.at("Size [bits]"), std::to_string(ibf.bit_size()));
    // Retired IBFs and deleted user bins only exist in an HIBF.
    EXPECT_FALSE(output.values.contains("Number of retired IBFs"));
    EXPECT_FALSE(output.values.contains("Number of deleted user bins"));
    EXPECT_TRUE(output.user_bins.empty());
    expect_text_ibfs(output, {expected});

    // The parameters are not affected by --detailed.
    cli_test_result const parameters_only = execute_app("raptor", "inspect", "--index", index);
    EXPECT_TRUE(result.out.starts_with(parameters_only.out));
}

TEST_P(inspect_ibf, json)
{
    auto const [number_of_repeated_bins, window_size] = GetParam();
    std::filesystem::path const index = ibf_path(number_of_repeated_bins, window_size);

    cli_test_result const result =
        execute_app("raptor", "inspect", "--index", index, "--format json", "--user-bins", "--detailed");
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);

    auto const loaded = load_index<ibf_index_t>(index);
    rapidjson::Document const json = parse_json(result.out);
    ASSERT_TRUE(json.IsObject());

    EXPECT_EQ(std::string{json["index"]["path"].GetString()}, index.string());
    EXPECT_EQ(json["index"]["version"].GetUint64(), 3u);
    EXPECT_EQ(std::string{json["index"]["type"].GetString()}, "IBF");
    EXPECT_EQ(json["index"]["parts"].GetUint64(), 1u);
    EXPECT_EQ(json["index"]["size_in_bytes"].GetUint64(), std::filesystem::file_size(index));

    EXPECT_EQ(json["minimiser"]["window_size"].GetUint64(), window_size);
    EXPECT_EQ(std::string{json["minimiser"]["shape"].GetString()}, "1111111111111111111");
    EXPECT_EQ(json["minimiser"]["shape_size"].GetUint64(), 19u);
    EXPECT_EQ(json["minimiser"]["shape_count"].GetUint64(), 19u);

    EXPECT_EQ(json["bloom_filter"]["false_positive_rate"].GetDouble(), 0.05);
    // These sections describe the HIBF.
    EXPECT_FALSE(json["bloom_filter"].HasMember("hash_functions"));
    EXPECT_FALSE(json.HasMember("layout"));
    EXPECT_FALSE(json.HasMember("update"));

    rapidjson::Value const & paths = json["user_bins"]["paths"];
    EXPECT_EQ(json["user_bins"]["count"].GetUint64(), loaded.bin_path().size());
    ASSERT_EQ(paths.Size(), loaded.bin_path().size());
    for (size_t user_bin_id = 0; user_bin_id < loaded.bin_path().size(); ++user_bin_id)
    {
        auto const & files = paths[static_cast<rapidjson::SizeType>(user_bin_id)];
        ASSERT_EQ(files.Size(), loaded.bin_path()[user_bin_id].size());
        for (size_t i = 0; i < files.Size(); ++i)
            EXPECT_EQ(std::string{files[static_cast<rapidjson::SizeType>(i)].GetString()},
                      loaded.bin_path()[user_bin_id][i]);
    }

    rapidjson::Value const & structure = json["structure"];
    EXPECT_EQ(structure["levels"].GetUint64(), 1u);
    EXPECT_EQ(structure["technical_bins"].GetUint64(), loaded.ibf().bin_count());
    EXPECT_EQ(structure["bits"].GetUint64(), loaded.ibf().bit_size());
    EXPECT_FALSE(structure.HasMember("retired_ibfs"));
    EXPECT_FALSE(structure.HasMember("deleted_user_bins"));
    expect_json_ibfs(structure["ibfs"], {expected_from(loaded.ibf(), 0u)});
}

INSTANTIATE_TEST_SUITE_P(inspect_ibf_suite,
                         inspect_ibf,
                         testing::Combine(testing::Values(0, 16, 32), testing::Values(19, 23)),
                         [](testing::TestParamInfo<inspect_ibf::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::max<int>(1, std::get<0>(info.param) * 4)) + "_bins_"
                                              + std::to_string(std::get<1>(info.param)) + "_window";
                             return name;
                         });

// The user bins are listed exactly as in the header of a `raptor search` result, including user bins with several
// files. The shipped indices have one file per user bin.
TEST_F(inspect_test, same_user_bins_as_search)
{
    std::string const bin1{data("bin1.fa").string()};
    std::string const bin2{data("bin2.fa").string()};
    std::string const bin3{data("bin3.fa").string()};
    std::string const bin4{data("bin4.fa").string()};
    {
        std::ofstream bin_file{"bins.txt"};
        bin_file << bin1 << ' ' << bin2 << '\n' << bin3 << '\n' << bin4 << ' ' << bin1 << ' ' << bin2 << '\n';
    }

    cli_test_result const build = execute_app("raptor",
                                              "build",
                                              "--quiet",
                                              "--kmer 19",
                                              "--window 19",
                                              "--input bins.txt",
                                              "--output raptor.index");
    ASSERT_EQ(build.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(build);

    cli_test_result const search = execute_app("raptor",
                                               "search",
                                               "--quiet",
                                               "--index raptor.index",
                                               "--query",
                                               data("query.fq"),
                                               "--output search.out");
    RAPTOR_ASSERT_ZERO_EXIT(search);

    std::vector<std::string> search_user_bins{};
    std::ifstream search_result{"search.out"};
    for (std::string line{}; std::getline(search_result, line) && line != "#QUERY_NAME\tUSER_BINS";)
        if (!line.starts_with("##"))
            search_user_bins.push_back(line);

    // The files of a user bin are separated by commas.
    EXPECT_EQ(
        search_user_bins,
        (std::vector<std::string>{"#0\t" + bin1 + ',' + bin2, "#1\t" + bin3, "#2\t" + bin4 + ',' + bin1 + ',' + bin2}));

    cli_test_result const text = execute_app("raptor", "inspect", "--index raptor.index", "--user-bins");
    EXPECT_EQ(text.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(text);
    EXPECT_EQ(parse_text(text.out).user_bins, search_user_bins);

    // JSON lists the files of a user bin as an array.
    cli_test_result const json =
        execute_app("raptor", "inspect", "--index raptor.index", "--user-bins", "--format json");
    EXPECT_EQ(json.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(json);

    rapidjson::Document const document = parse_json(json.out);
    ASSERT_TRUE(document.IsObject());
    rapidjson::Value const & paths = document["user_bins"]["paths"];
    std::vector<std::vector<std::string>> json_paths{};
    for (auto const & files : paths.GetArray())
    {
        json_paths.emplace_back();
        for (auto const & file : files.GetArray())
            json_paths.back().emplace_back(file.GetString());
    }
    EXPECT_EQ(json_paths, (std::vector<std::vector<std::string>>{{bin1, bin2}, {bin3}, {bin4, bin1, bin2}}));
}

struct inspect_ibf_partitioned : public inspect_test
{
    static constexpr size_t parts{4u};

    void SetUp() override
    {
        inspect_test::SetUp();

        std::ofstream bin_file{"bins.txt"};
        for (auto && file_path : get_repeated_bins(16u))
            bin_file << file_path << '\n';
        bin_file.close();

        cli_test_result const build = execute_app("raptor",
                                                  "build",
                                                  "--quiet",
                                                  "--kmer 19",
                                                  "--window 23",
                                                  "--parts",
                                                  std::to_string(parts),
                                                  "--input bins.txt",
                                                  "--output raptor.index");
        ASSERT_EQ(build.err, std::string{});
        RAPTOR_ASSERT_ZERO_EXIT(build);
    }

    static std::string part(size_t const id)
    {
        return "raptor.index_" + std::to_string(id);
    }
};

TEST_F(inspect_ibf_partitioned, detailed)
{
    cli_test_result const result = execute_app("raptor", "inspect", "--index raptor.index", "--detailed");
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);

    size_t size_in_bytes{};
    size_t bits{};
    size_t technical_bins{};
    std::vector<expected_ibf> expected{};
    for (size_t id = 0; id < parts; ++id)
    {
        auto const loaded = load_index<ibf_index_t>(part(id));
        size_in_bytes += std::filesystem::file_size(part(id));
        bits += loaded.ibf().bit_size();
        technical_bins += loaded.ibf().bin_count();
        expected.push_back(expected_from(loaded.ibf(), id));
    }

    text_output const output = parse_text(result.out);
    EXPECT_EQ(output.values.at("Index"), "\"raptor.index\"");
    EXPECT_EQ(output.values.at("Index parts"), std::to_string(parts));
    EXPECT_EQ(output.values.at("Index size [bytes]"), std::to_string(size_in_bytes));
    EXPECT_EQ(output.values.at("Window size"), "23");
    EXPECT_EQ(output.values.at("Number of user bins"), "64");
    EXPECT_EQ(output.values.at("Levels"), "1");
    EXPECT_EQ(output.values.at("Technical bins"), std::to_string(technical_bins));
    EXPECT_EQ(output.values.at("Size [bits]"), std::to_string(bits));
    expect_text_ibfs(output, expected);
}

TEST_F(inspect_ibf_partitioned, part_instead_of_index)
{
    cli_test_result const result = execute_app("raptor", "inspect", "--index", part(0u));
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err,
              std::string{"[Error] The index \"raptor.index_0\" is a part of a partitioned index with 4 parts. Pass "
                          "the path without the suffix _<part>.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(inspect_ibf_partitioned, missing_part)
{
    std::filesystem::remove(part(parts - 1u));

    cli_test_result const result = execute_app("raptor", "inspect", "--index raptor.index");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, "[Error] The file \"" + part(parts - 1u) + "\" does not exist!\n");
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

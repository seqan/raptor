// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "inspect_test.hpp"

struct inspect_hibf : public inspect_test, public testing::WithParamInterface<std::tuple<size_t, size_t>>
{};

// The configuration of the shipped HIBF indices is the one in test/data/hibf/<N>bins.layout.
TEST_P(inspect_hibf, text)
{
    auto const [number_of_repeated_bins, window_size] = GetParam();
    std::filesystem::path const index = ibf_path(number_of_repeated_bins, window_size, is_hibf::yes);
    std::string const number_of_user_bins = std::to_string(std::max<size_t>(1u, number_of_repeated_bins * 4u));
    std::string const number_of_ibfs = std::to_string(load_index<hibf_index_t>(index).ibf().ibf_vector.size());

    cli_test_result const result = execute_app("raptor", "inspect", "--index", index);
    std::string const expected = "### Index\n"
                                 "## Index = \""
                               + index.string()
                               + "\"\n"
                                 "## Index version = 3\n"
                                 "## Index type = HIBF\n"
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
                                 "## Relaxed false positive rate = 0.05\n"
                                 "## Hash functions = 2\n"
                                 "### Layout parameters\n"
                                 "## Number of user bins in layout = "
                               + number_of_user_bins
                               + "\n"
                                 "## tmax = 64\n"
                                 "## Empty bin fraction = 0\n"
                                 "## Sketch bits = 12\n"
                                 "## alpha = 1.2\n"
                                 "## Max rearrangement ratio = 0.5\n"
                                 "## Disable estimate union = false\n"
                                 "## Disable rearrangement = false\n"
                                 "### Update state\n"
                                 "## Track occupancy = false\n"
                                 "## Number of IBFs = "
                               + number_of_ibfs
                               + "\n"
                                 "## Original number of IBFs = "
                               + number_of_ibfs
                               + "\n"
                                 "## Resized IBFs = 0\n"
                                 "### User bins\n"
                                 "## Number of user bins = "
                               + number_of_user_bins + "\n";
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);
}

INSTANTIATE_TEST_SUITE_P(inspect_hibf_suite,
                         inspect_hibf,
                         testing::Combine(testing::Values(0, 16, 32), testing::Values(19, 23)),
                         [](testing::TestParamInfo<inspect_hibf::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::max<int>(1, std::get<0>(info.param) * 4)) + "_bins_"
                                              + std::to_string(std::get<1>(info.param)) + "_window";
                             return name;
                         });

TEST_F(inspect_test, three_levels_text)
{
    std::filesystem::path const index = data("three_levels.hibf");
    cli_test_result const result = execute_app("raptor", "inspect", "--index", index, "--detailed", "--user-bins");
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);

    hibf_index_t const loaded = load_index<hibf_index_t>(index);
    std::vector<expected_ibf> const expected = expected_ibfs(loaded);

    // The layout has 49 IBFs on three levels, each with 8 technical bins.
    ASSERT_EQ(expected.size(), 49u);
    EXPECT_EQ(std::ranges::max(expected, {}, &expected_ibf::level).level, 2u);

    text_output const output = parse_text(result.out);
    EXPECT_EQ(output.sections,
              (std::vector<std::string>{"Index",
                                        "Minimiser parameters",
                                        "Bloom filter parameters",
                                        "Layout parameters",
                                        "Update state",
                                        "User bins",
                                        "Structure"}));
    EXPECT_EQ(output.values.at("Relaxed false positive rate"), "0.3");
    EXPECT_EQ(output.values.at("Levels"), "3");
    EXPECT_EQ(output.values.at("Technical bins"), std::to_string(49u * 8u));
    EXPECT_EQ(output.values.at("Number of retired IBFs"), "0");
    EXPECT_EQ(output.values.at("Number of deleted user bins"), "0");
    // Only listed if there are any.
    EXPECT_FALSE(output.values.contains("Retired IBFs"));
    EXPECT_FALSE(output.values.contains("Deleted user bins"));
    EXPECT_EQ(output.user_bins, user_bin_lines(loaded.bin_path()));
    expect_text_ibfs(output, expected);

    // As in the JSON output, the files belong to the user bins section and precede the structure.
    EXPECT_NE(result.out.find("## Number of user bins = 128\n#0\t"), std::string::npos);
    EXPECT_LT(result.out.find("\n#127\t"), result.out.find("### Structure\n"));
}

TEST_F(inspect_test, three_levels_json)
{
    std::filesystem::path const index = data("three_levels.hibf");
    cli_test_result const result = execute_app("raptor", "inspect", "--index", index, "--detailed", "--format json");
    EXPECT_EQ(result.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result);

    hibf_index_t const loaded = load_index<hibf_index_t>(index);
    auto const & config = loaded.config();
    rapidjson::Document const json = parse_json(result.out);
    ASSERT_TRUE(json.IsObject());

    EXPECT_EQ(std::string{json["index"]["type"].GetString()}, "HIBF");
    EXPECT_EQ(json["bloom_filter"]["false_positive_rate"].GetDouble(), loaded.fpr());
    EXPECT_EQ(json["bloom_filter"]["relaxed_false_positive_rate"].GetDouble(), config.relaxed_fpr);
    EXPECT_EQ(json["bloom_filter"]["hash_functions"].GetUint64(), config.number_of_hash_functions);

    rapidjson::Value const & layout = json["layout"];
    EXPECT_EQ(layout["number_of_user_bins"].GetUint64(), config.number_of_user_bins);
    EXPECT_EQ(layout["tmax"].GetUint64(), config.tmax);
    EXPECT_EQ(layout["empty_bin_fraction"].GetDouble(), config.empty_bin_fraction);
    EXPECT_EQ(layout["sketch_bits"].GetUint64(), config.sketch_bits);
    EXPECT_EQ(layout["alpha"].GetDouble(), config.alpha);
    EXPECT_EQ(layout["max_rearrangement_ratio"].GetDouble(), config.max_rearrangement_ratio);
    EXPECT_EQ(layout["disable_estimate_union"].GetBool(), config.disable_estimate_union);
    EXPECT_EQ(layout["disable_rearrangement"].GetBool(), config.disable_rearrangement);

    rapidjson::Value const & update = json["update"];
    EXPECT_EQ(update["track_occupancy"].GetBool(), config.track_occupancy);
    EXPECT_EQ(update["number_of_ibfs"].GetUint64(), loaded.ibf().ibf_vector.size());
    EXPECT_EQ(update["original_number_of_ibfs"].GetUint64(), loaded.original_number_of_ibfs());
    EXPECT_EQ(update["resized_ibfs"].GetUint64(), loaded.number_of_resized_ibfs());

    EXPECT_EQ(json["user_bins"]["count"].GetUint64(), 128u);
    EXPECT_FALSE(json["user_bins"].HasMember("paths"));

    rapidjson::Value const & structure = json["structure"];
    EXPECT_EQ(structure["levels"].GetUint64(), 3u);
    EXPECT_EQ(structure["technical_bins"].GetUint64(), 49u * 8u);
    EXPECT_TRUE(json_ids(structure["retired_ibfs"]).empty());
    EXPECT_TRUE(json_ids(structure["deleted_user_bins"]).empty());
    expect_json_ibfs(structure["ibfs"], expected_ibfs(loaded));
}

/*!\brief An index that has been updated: grown by insertions, and shrunk by deletions that retire an IBF.
 * \details
 * The shipped indices have no empty bins and cannot be updated, so this builds its own. The assertions on the
 * loaded index make sure that the workload still produces resized IBFs, retired IBFs and deleted user bins; the
 * output of `raptor inspect` is compared against the loaded index.
 */
struct inspect_dynamic : public inspect_test
{
    std::vector<std::string> paths{};

    void SetUp() override
    {
        inspect_test::SetUp();
        paths = generate_bins("bins", 140u);

        std::vector<std::string> const initial(paths.begin(), paths.begin() + 120);
        std::vector<std::string> const inserted(paths.begin() + 120, paths.end());
        ASSERT_NO_FATAL_FAILURE(build_index(write_bin_file("initial.txt", initial), "built.index"));

        cli_test_result const insertion = execute_app("raptor",
                                                      "update",
                                                      "insert",
                                                      "--index built.index",
                                                      "--output inserted.index",
                                                      "--insert",
                                                      write_bin_file("inserted.txt", inserted));
        ASSERT_EQ(insertion.err, std::string{});
        RAPTOR_ASSERT_ZERO_EXIT(insertion);
    }
};

TEST_F(inspect_dynamic, update_state)
{
    hibf_index_t const built = load_index<hibf_index_t>("built.index");
    hibf_index_t const inserted = load_index<hibf_index_t>("inserted.index");
    ASSERT_GT(inserted.number_of_resized_ibfs(), 0u) << "The workload no longer resizes an IBF.";

    for (auto const & [file, loaded] : {std::pair{"built.index", &built}, std::pair{"inserted.index", &inserted}})
    {
        cli_test_result const result = execute_app("raptor", "inspect", "--index", file);
        EXPECT_EQ(result.err, std::string{});
        RAPTOR_ASSERT_ZERO_EXIT(result);

        text_output const output = parse_text(result.out);
        // The layout is computed with `--empty-bin-fraction 0.1`. seqan::hibf::config::validate_and_set_defaults
        // rounds it to whole technical bins of `tmax`: 57 of 64 bins hold data, i.e. 7/64 are empty.
        EXPECT_EQ(output.values.at("Empty bin fraction"), "0.109375") << file;
        EXPECT_EQ(output.values.at("Track occupancy"), "true") << file;
        // Read via load_parameters, compared against the whole index.
        EXPECT_EQ(output.values.at("Number of IBFs"), std::to_string(loaded->ibf().ibf_vector.size())) << file;
        EXPECT_EQ(output.values.at("Original number of IBFs"), std::to_string(loaded->original_number_of_ibfs()))
            << file;
        EXPECT_EQ(output.values.at("Resized IBFs"), std::to_string(loaded->number_of_resized_ibfs())) << file;
        EXPECT_EQ(output.values.at("Number of user bins"), std::to_string(loaded->bin_path().size())) << file;
    }

    // Inserting does not change the layout parameters.
    cli_test_result const result = execute_app("raptor", "inspect", "--index inserted.index");
    text_output const output = parse_text(result.out);
    EXPECT_EQ(output.values.at("Number of user bins in layout"), "120");
    EXPECT_EQ(output.values.at("Number of user bins"), "140");
}

TEST_F(inspect_dynamic, deleted_user_bins_and_retired_ibfs)
{
    // Deleting every user bin of a lower-level IBF without merged bins retires it. User bin 0 is deleted as well;
    // unless it belongs to that IBF, this only frees its technical bins.
    std::set<size_t> user_bins_to_delete{0u};
    {
        hibf_index_t const inserted = load_index<hibf_index_t>("inserted.index");
        auto const & ibf_bin_to_user_bin_id = inserted.ibf().ibf_bin_to_user_bin_id;
        bool found{false};
        for (size_t ibf_idx = 1u; ibf_idx < ibf_bin_to_user_bin_id.size() && !found; ++ibf_idx)
        {
            auto const & user_bin_ids = ibf_bin_to_user_bin_id[ibf_idx];
            if (std::ranges::find(user_bin_ids, seqan::hibf::bin_kind::merged) == user_bin_ids.end())
            {
                for (uint64_t const user_bin_id : user_bin_ids)
                    if (user_bin_id < inserted.bin_path().size()) // Not seqan::hibf::bin_kind::deleted.
                        user_bins_to_delete.insert(user_bin_id);
                found = true;
            }
        }
        ASSERT_TRUE(found) << "There is no lower-level IBF without merged bins.";
    }

    std::string delete_options{};
    for (size_t const user_bin_id : user_bins_to_delete)
        delete_options += "--delete " + std::to_string(user_bin_id) + ' ';

    cli_test_result const deletion =
        execute_app("raptor", "update", "delete", "--index inserted.index", "--output deleted.index", delete_options);
    ASSERT_EQ(deletion.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(deletion);

    hibf_index_t const deleted = load_index<hibf_index_t>("deleted.index");
    std::vector<size_t> const retired_ibfs = expected_retired_ibfs(deleted);
    std::vector<size_t> const deleted_user_bins = expected_deleted_user_bins(deleted);
    std::vector<expected_ibf> const ibfs = expected_ibfs(deleted);
    ASSERT_FALSE(retired_ibfs.empty()) << "The workload no longer retires an IBF.";
    ASSERT_EQ(deleted_user_bins, std::vector<size_t>(user_bins_to_delete.begin(), user_bins_to_delete.end()));

    // Text
    {
        cli_test_result const result = execute_app("raptor", "inspect", "--index deleted.index", "--detailed");
        EXPECT_EQ(result.err, std::string{});
        RAPTOR_ASSERT_ZERO_EXIT(result);

        text_output const output = parse_text(result.out);
        // Deleting keeps the IDs: The user bins are still counted, and retired IBFs still occupy their slot.
        EXPECT_EQ(output.values.at("Number of user bins"), "140");
        EXPECT_EQ(output.values.at("Number of IBFs"), std::to_string(deleted.ibf().ibf_vector.size()));
        EXPECT_EQ(output.values.at("Number of retired IBFs"), std::to_string(retired_ibfs.size()));
        EXPECT_EQ(output.values.at("Retired IBFs"), comma_separated(retired_ibfs));
        EXPECT_EQ(output.values.at("Number of deleted user bins"), std::to_string(deleted_user_bins.size()));
        EXPECT_EQ(output.values.at("Deleted user bins"), comma_separated(deleted_user_bins));
        expect_text_ibfs(output, ibfs);
    }

    // JSON
    {
        cli_test_result const result =
            execute_app("raptor", "inspect", "--index deleted.index", "--detailed", "--format json");
        EXPECT_EQ(result.err, std::string{});
        RAPTOR_ASSERT_ZERO_EXIT(result);

        rapidjson::Document const json = parse_json(result.out);
        ASSERT_TRUE(json.IsObject());
        EXPECT_EQ(json_ids(json["structure"]["retired_ibfs"]), retired_ibfs);
        EXPECT_EQ(json_ids(json["structure"]["deleted_user_bins"]), deleted_user_bins);
        expect_json_ibfs(json["structure"]["ibfs"], ibfs);
    }
}

// --detailed checks that the structure of the HIBF is consistent. Without it, the structure is not read at all.
TEST_F(inspect_test, inconsistent_parent)
{
    hibf_index_t index = load_index<hibf_index_t>(data("three_levels.hibf"));
    index.ibf().prev_ibf_id[1].ibf_idx = 1u; // A cycle.
    store_index("broken.index", index);

    cli_test_result const parameters_only = execute_app("raptor", "inspect", "--index broken.index");
    RAPTOR_ASSERT_ZERO_EXIT(parameters_only);

    cli_test_result const result = execute_app("raptor", "inspect", "--index broken.index", "--detailed");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[Error] Cannot read index: IBF 1 is not connected to the top-level IBF.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(inspect_test, inconsistent_user_bin)
{
    hibf_index_t index = load_index<hibf_index_t>(data("three_levels.hibf"));
    ASSERT_NE(index.ibf().ibf_bin_to_user_bin_id[2][0], seqan::hibf::bin_kind::merged);
    index.ibf().ibf_bin_to_user_bin_id[2][0] = 1000u;
    store_index("broken.index", index);

    cli_test_result const result = execute_app("raptor", "inspect", "--index broken.index", "--detailed");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err,
              std::string{"[Error] Cannot read index: IBF 2 refers to user bin 1000, but the index has only 128 user "
                          "bins.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

TEST_F(inspect_test, inconsistent_sizes)
{
    hibf_index_t index = load_index<hibf_index_t>(data("three_levels.hibf"));
    index.ibf().ibf_bin_to_user_bin_id.pop_back();
    store_index("broken.index", index);

    cli_test_result const result = execute_app("raptor", "inspect", "--index broken.index", "--detailed");
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err,
              std::string{"[Error] Cannot read index: The HIBF does not store the same number of entries for every "
                          "IBF.\n"});
    RAPTOR_ASSERT_FAIL_EXIT(result);
}

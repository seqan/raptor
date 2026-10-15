// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <cmath>
#include <cstdint>
#include <limits>
#include <sstream>
#include <string>

#include <cereal/archives/json.hpp> // Configures the rapidjson bundled with cereal.

#include <raptor/inspect/json_writer.hpp>

// Strict RFC 8259: Unlike cereal's default flags, this rejects NaN and infinity.
static bool is_valid_json(std::string const & text)
{
    rapidjson::Document document{};
    document.Parse<rapidjson::kParseFullPrecisionFlag | rapidjson::kParseValidateEncodingFlag>(text.data(),
                                                                                               text.size());
    return !document.HasParseError();
}

TEST(json_quote, plain)
{
    EXPECT_EQ(raptor::json_quote(""), "\"\"");
    EXPECT_EQ(raptor::json_quote("/path/to/bin1.fa"), "\"/path/to/bin1.fa\"");
}

TEST(json_quote, escapes)
{
    EXPECT_EQ(raptor::json_quote("a\"b"), R"("a\"b")");
    EXPECT_EQ(raptor::json_quote("a\\b"), R"("a\\b")");
    EXPECT_EQ(raptor::json_quote("\b\f\n\r\t"), R"("\b\f\n\r\t")");
}

TEST(json_quote, control_characters)
{
    using namespace std::string_literals;
    EXPECT_EQ(raptor::json_quote("\0"s), R"("\u0000")");
    EXPECT_EQ(raptor::json_quote("\x01"), R"("\u0001")");
    EXPECT_EQ(raptor::json_quote("\x1f"), R"("\u001f")");
    // Neither 0x20 nor DEL (0x7f) need to be escaped.
    EXPECT_EQ(raptor::json_quote(" \x7f"), "\" \x7f\"");
}

// Bytes outside ASCII are passed through, i.e. UTF-8 stays UTF-8.
TEST(json_quote, utf8)
{
    EXPECT_EQ(raptor::json_quote("Gänsefüßchen"), "\"Gänsefüßchen\"");
}

TEST(shortest_double, values)
{
    EXPECT_EQ(raptor::shortest_double(0.05), "0.05");
    EXPECT_EQ(raptor::shortest_double(1.2), "1.2");
    EXPECT_EQ(raptor::shortest_double(0.0), "0");
    EXPECT_EQ(raptor::shortest_double(1.0), "1");
    EXPECT_EQ(raptor::shortest_double(1e-05), "1e-05");
    EXPECT_EQ(raptor::shortest_double(0.1 + 0.2), "0.30000000000000004");
    EXPECT_EQ(raptor::shortest_double(std::numeric_limits<double>::lowest()), "-1.7976931348623157e+308");
}

TEST(json_writer, empty_containers)
{
    std::ostringstream stream{};
    raptor::json_writer json{stream};
    json.begin_object();
    json.key("object").begin_object().end_object();
    json.key("array").begin_array().end_array();
    json.key("inline_array").begin_array(true).end_array();
    json.end_object();

    EXPECT_EQ(stream.str(),
              "{\n"
              "  \"object\": {},\n"
              "  \"array\": [],\n"
              "  \"inline_array\": []\n"
              "}\n");
    EXPECT_TRUE(is_valid_json(stream.str()));
}

TEST(json_writer, values)
{
    std::ostringstream stream{};
    raptor::json_writer json{stream};
    json.begin_object();
    json.key("string").string("text");
    json.key("true").boolean(true);
    json.key("false").boolean(false);
    json.key("null").null();
    // An `uint8_t` must be written as a number, not as a character.
    json.key("uint8_t").number(uint8_t{19u});
    json.key("uint64_t").number(std::numeric_limits<uint64_t>::max());
    json.key("double").number(0.05);
    json.key("float").number(0.5f);
    json.key("nan").number(std::numeric_limits<double>::quiet_NaN());
    json.key("infinity").number(-std::numeric_limits<double>::infinity());
    json.end_object();

    EXPECT_EQ(stream.str(),
              "{\n"
              "  \"string\": \"text\",\n"
              "  \"true\": true,\n"
              "  \"false\": false,\n"
              "  \"null\": null,\n"
              "  \"uint8_t\": 19,\n"
              "  \"uint64_t\": 18446744073709551615,\n"
              "  \"double\": 0.05,\n"
              "  \"float\": 0.5,\n"
              "  \"nan\": null,\n"
              "  \"infinity\": null\n"
              "}\n");
    EXPECT_TRUE(is_valid_json(stream.str()));
}

TEST(json_writer, nesting)
{
    std::ostringstream stream{};
    raptor::json_writer json{stream};
    json.begin_object();
    json.key("outer").begin_object();
    json.key("inner").begin_array();
    json.number(1u);
    json.begin_object().key("a").number(2u).end_object();
    json.end_array();
    json.end_object();
    json.key("key \"with\" quotes").string("value\n");
    json.end_object();

    EXPECT_EQ(stream.str(),
              "{\n"
              "  \"outer\": {\n"
              "    \"inner\": [\n"
              "      1,\n"
              "      {\n"
              "        \"a\": 2\n"
              "      }\n"
              "    ]\n"
              "  },\n"
              "  \"key \\\"with\\\" quotes\": \"value\\n\"\n"
              "}\n");
    EXPECT_TRUE(is_valid_json(stream.str()));
}

// Everything nested in an inline container is written on the same line.
TEST(json_writer, inline_elements)
{
    std::ostringstream stream{};
    raptor::json_writer json{stream};
    json.begin_array();
    json.begin_array(true).string("a").string("b").end_array();
    json.begin_object(true);
    json.key("id").number(0u);
    json.key("children").begin_array().number(1u).number(2u).end_array();
    json.key("nested").begin_object().key("x").null().end_object();
    json.end_object();
    json.end_array();

    EXPECT_EQ(stream.str(),
              "[\n"
              "  [\"a\", \"b\"],\n"
              "  {\"id\": 0, \"children\": [1, 2], \"nested\": {\"x\": null}}\n"
              "]\n");
    EXPECT_TRUE(is_valid_json(stream.str()));
}

TEST(json_writer, top_level_array)
{
    std::ostringstream stream{};
    raptor::json_writer json{stream};
    json.begin_array(true).number(1u).end_array();

    EXPECT_EQ(stream.str(), "[1]\n");
    EXPECT_TRUE(is_valid_json(stream.str()));
}

// Confirms that the oracle rejects what the writer must never produce.
TEST(json_writer, oracle)
{
    EXPECT_FALSE(is_valid_json("{\"nan\": NaN}"));
    EXPECT_FALSE(is_valid_json("[1, 2,]"));
    EXPECT_FALSE(is_valid_json("[\"\x01\"]"));
    EXPECT_FALSE(is_valid_json("[\"\xff\"]"));
}

#include "utils.hpp"

#include "selections/SelectionParser.hpp"
#include "selections/boolean.hpp"
#include "selections/properties.hpp"
#include <molpp/internal/MolData.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace mol::internal;
using namespace testing;

//! Test fixture for SelectionParser selection nodes.
class SelectionParserTest : public testing::Test
{
public:
    SelectionParserTest()
    : parser{default_parser()}
    {
    }

    SelectionParser parser;
};

// Constructing with an invalid grammar should throw an exception.
TEST_F(SelectionParserTest, invalid_grammar)
{
    EXPECT_THROW(SelectionParser{""}, mol::Error);
}

// Parsing an empty string should throw an exception.
TEST_F(SelectionParserTest, empty_expression)
{
    EXPECT_THROW(parser.parse(""), mol::Error);
}

// Parsing an invalid expression should throw an exception.
TEST_F(SelectionParserTest, invalid_expression)
{
    EXPECT_THROW(parser.parse("invalid"), mol::Error);
}

// Parsing an invalid expression should throw an exception with an error message.
TEST_F(SelectionParserTest, invalid_expression_message)
{
    std::string message;

    try
    {
        parser.parse("invalid");
    }
    catch (mol::Error const& error)
    {
        message = error.what();
    }

    EXPECT_THAT(message, HasSubstr("invalid"));
    EXPECT_THAT(message, HasSubstr("^"));
}

// Parsing a valid expression should return a non-null selection node.
TEST_F(SelectionParserTest, valid_expression)
{
    std::shared_ptr<SelectionNode> const node = parser.parse("resid 1");

    EXPECT_THAT(node, NotNull());
}

// Parsing an expression containing "or" should return a valid selection node tree.
TEST_F(SelectionParserTest, or_expression)
{
    std::shared_ptr<SelectionNode> const tree = parser.parse("resid 1 or resid 2 or resid 3");

    ASSERT_TRUE(std::dynamic_pointer_cast<OrSelection>(tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<OrSelection>(tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(tree->left->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(tree->left->right));
}

// Parsing an expression containing "and" should return a valid selection node tree.
TEST_F(SelectionParserTest, and_expression)
{
    std::shared_ptr<SelectionNode> const tree = parser.parse("resid 1 and resid 2 and resid 3");

    ASSERT_TRUE(std::dynamic_pointer_cast<AndSelection>(tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<AndSelection>(tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(tree->left->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(tree->left->right));
}

// Parsing an expression mixding "and" and "or" should return a valid selection node tree.
TEST_F(SelectionParserTest, and_or_expression)
{
    std::shared_ptr<SelectionNode> const tree = parser.parse("resid 1 or resid 2 and resid 3");

    ASSERT_TRUE(std::dynamic_pointer_cast<AndSelection>(tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<OrSelection>(tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(tree->left->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(tree->left->right));
}

// Parsing an expression containing "not" should return a valid selection node tree.
TEST_F(SelectionParserTest, not_expression)
{
    std::shared_ptr<SelectionNode> const tree = parser.parse("not resid 1 or resid 2");

    ASSERT_TRUE(std::dynamic_pointer_cast<OrSelection>(tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<NotSelection>(tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(tree->left->left));
}

// Parsing an expression containing "not" between boolean operators should return a valid selection node tree.
TEST_F(SelectionParserTest, not_between_boolean_expression)
{
    std::shared_ptr<SelectionNode> const tree = parser.parse("resid 1 and not resid 2");

    ASSERT_TRUE(std::dynamic_pointer_cast<AndSelection>(tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<NotSelection>(tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(tree->right->left));
}

// Parsing an expression containing parentheses should return a valid selection node tree.
TEST_F(SelectionParserTest, parentheses_expression)
{
    std::shared_ptr<SelectionNode> const tree = parser.parse("resid 1 or not (resid 2 or resid 3)");

    ASSERT_TRUE(std::dynamic_pointer_cast<OrSelection>(tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<NotSelection>(tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<OrSelection>(tree->right->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(tree->right->left->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(tree->right->left->right));
}

// Parsing an expression containing nested parentheses should return a valid selection node tree.
TEST_F(SelectionParserTest, nested_parentheses_expression)
{
    std::shared_ptr<SelectionNode> const tree = parser.parse("resid 1 or (resid 2 or (resid 3 and resid 4))");

    ASSERT_TRUE(std::dynamic_pointer_cast<OrSelection>(tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<OrSelection>(tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(tree->right->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<AndSelection>(tree->right->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(tree->right->right->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(tree->right->right->right));
}

// Parsing an expression containing "all" should return a valid selection node tree.
TEST_F(SelectionParserTest, all_expression)
{
    std::shared_ptr<SelectionNode> const tree = parser.parse("all");

    ASSERT_TRUE(std::dynamic_pointer_cast<AllSelection>(tree));
}

// Parsing an expression containing "all" and other operators should return a valid selection node tree.
TEST_F(SelectionParserTest, all_and_expression)
{
    std::shared_ptr<SelectionNode> const tree = parser.parse("all and resid 1");

    ASSERT_TRUE(std::dynamic_pointer_cast<AndSelection>(tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<AllSelection>(tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(tree->right));
}

// Parsing a resid expression with one number should return a valid selection node tree.
TEST_F(SelectionParserTest, resid_number)
{
    std::shared_ptr<SelectionNode> const tree = parser.parse("resid 1");

    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(tree));
}

// Parsing a resid expression with a range should return a valid selection node tree.
TEST_F(SelectionParserTest, resid_range)
{
    std::shared_ptr<SelectionNode> const tree = parser.parse("resid 1:10");

    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(tree));
}

// Parsing a resid expression with a range and a single number should return a valid selection node tree.
TEST_F(SelectionParserTest, resid_range_number)
{
    std::shared_ptr<SelectionNode> const tree = parser.parse("resid 1:10 20");

    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(tree));
}

// Parsing a resid expression with a number, a range and a number should return a valid selection node tree.
TEST_F(SelectionParserTest, resid_number_range_number)
{
    std::shared_ptr<SelectionNode> const tree = parser.parse("resid 1 2:10 20");

    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(tree));
}

// Parsing a resid expression with multiple numbers and ranges should return a valid selection node tree.
TEST_F(SelectionParserTest, resid_multiple_numbers_ranges)
{
    std::shared_ptr<SelectionNode> const tree = parser.parse("resid 1 2:10 20 30:40");

    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(tree));
}

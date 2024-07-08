#include <molpp/AtomSelector.hpp>
#include "selections/boolean.hpp"
#include "selections/properties.hpp"
#include "selections/SelectionStack.hpp"
#include "selections/SelectionParser.hpp"
#include <molpp/Error.hpp>

#include "utils.hpp"
#include "files.hpp"

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <unordered_set>

using namespace mol;
using namespace mol::internal;
using namespace testing;

std::shared_ptr<std::unordered_set<index_t>> evaluate_sel_tree(std::shared_ptr<SelectionNode> tree, MolData const& data)
{
    SelectionIndices flags;
    for (index_t atom_idx = 0; atom_idx < data.size<Atom>(); atom_idx++)
    {
        flags.available->insert(atom_idx);
    }

    SelectionStack stack;
    stack.push_node(tree);
    stack.push_indices(flags);
    while (!stack.empty_nodes())
    {
        std::shared_ptr<SelectionNode> node = stack.pop_node();
        node->evaluate(stack, data, 0);
    }
    return flags.selected;
}

TEST(Selection, DataNodes) {
    SelNumber float_num("1.0001");
    EXPECT_EQ(float_num, 1.0001);
    EXPECT_EQ(float_num, 1.0002);
    EXPECT_NE(float_num, 1.0003);

    SelNumber int_num("1");
    EXPECT_EQ(int_num, 1);
    EXPECT_NE(int_num, 2);

    SelNumberRange range("3", "5");
    // Integer
    EXPECT_FALSE(range.has(-1));
    EXPECT_FALSE(range.has(2));
    EXPECT_TRUE(range.has(3));
    EXPECT_TRUE(range.has(4));
    EXPECT_TRUE(range.has(5));
    EXPECT_FALSE(range.has(6));
    // Floating point
    EXPECT_FALSE(range.has(-1.0));
    EXPECT_FALSE(range.has(2.0));
    EXPECT_TRUE(range.has(3.0));
    EXPECT_TRUE(range.has(4.0));
    EXPECT_TRUE(range.has(5.0));
    EXPECT_FALSE(range.has(6.0));

    SelNumberRange consecutive_range("3", "4");
    EXPECT_FALSE(consecutive_range.has(2));
    EXPECT_TRUE(consecutive_range.has(3));
    EXPECT_TRUE(consecutive_range.has(4));
    EXPECT_FALSE(consecutive_range.has(5));
}

TEST(Selection, SelectionParser) {
    EXPECT_THROW(SelectionParser("not valid"), Error);
    EXPECT_THROW(SEL_PARSER.parse("not valid"), Error);
    EXPECT_TRUE(SEL_PARSER.parse("resid 1"));
}

TEST(Selection, BooleanParsing) {
    // Mock MolData
    MolData data = create_moldata(5, 2, 1, 1, 1);
    std::shared_ptr<std::unordered_set<index_t>> selected;

    auto sel_tree = SEL_PARSER.parse("resid 1 or resid 2 or resid 3");
    ASSERT_TRUE(std::dynamic_pointer_cast<OrSelection>(sel_tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<OrSelection>(sel_tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->left->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->left->right));
    selected = evaluate_sel_tree(sel_tree, data);
    EXPECT_THAT(*selected, UnorderedElementsAre(2, 3, 4, 5, 6, 7));

    sel_tree = SEL_PARSER.parse("resid 1 and resid 2 and resid 3");
    ASSERT_TRUE(std::dynamic_pointer_cast<AndSelection>(sel_tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<AndSelection>(sel_tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->left->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->left->right));
    selected = evaluate_sel_tree(sel_tree, data);
    EXPECT_THAT(*selected, UnorderedElementsAre());

    sel_tree = SEL_PARSER.parse("resid 1 or resid 2 and resid 3");
    ASSERT_TRUE(std::dynamic_pointer_cast<AndSelection>(sel_tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<OrSelection>(sel_tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->left->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->left->right));
    selected = evaluate_sel_tree(sel_tree, data);
    EXPECT_THAT(*selected, UnorderedElementsAre());

    sel_tree = SEL_PARSER.parse("(resid 1 or resid 2) and resid 3");
    ASSERT_TRUE(std::dynamic_pointer_cast<AndSelection>(sel_tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<OrSelection>(sel_tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->left->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->left->right));
    selected = evaluate_sel_tree(sel_tree, data);
    EXPECT_THAT(*selected, UnorderedElementsAre());

    sel_tree = SEL_PARSER.parse("resid 1 or (resid 2 and resid 3)");
    ASSERT_TRUE(std::dynamic_pointer_cast<OrSelection>(sel_tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<AndSelection>(sel_tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->right->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->right->right));
    selected = evaluate_sel_tree(sel_tree, data);
    EXPECT_THAT(*selected, UnorderedElementsAre(2, 3));

    sel_tree = SEL_PARSER.parse("resid 1 and resid 2 or resid 3");
    ASSERT_TRUE(std::dynamic_pointer_cast<OrSelection>(sel_tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<AndSelection>(sel_tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->left->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->left->right));
    selected = evaluate_sel_tree(sel_tree, data);
    EXPECT_THAT(*selected, UnorderedElementsAre(6, 7));

    sel_tree = SEL_PARSER.parse("(resid 1 and resid 2) or resid 3");
    ASSERT_TRUE(std::dynamic_pointer_cast<OrSelection>(sel_tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<AndSelection>(sel_tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->left->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->left->right));
    selected = evaluate_sel_tree(sel_tree, data);
    EXPECT_THAT(*selected, UnorderedElementsAre(6, 7));

    sel_tree = SEL_PARSER.parse("resid 1 and (resid 2 or resid 3)");
    ASSERT_TRUE(std::dynamic_pointer_cast<AndSelection>(sel_tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<OrSelection>(sel_tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->right->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->right->right));
    selected = evaluate_sel_tree(sel_tree, data);
    EXPECT_THAT(*selected, UnorderedElementsAre());

    sel_tree = SEL_PARSER.parse("not resid 1 or resid 2");
    ASSERT_TRUE(std::dynamic_pointer_cast<OrSelection>(sel_tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<NotSelection>(sel_tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->left->left));
    selected = evaluate_sel_tree(sel_tree, data);
    EXPECT_THAT(*selected, UnorderedElementsAre(0, 1, 4, 5, 6, 7, 8, 9));

    sel_tree = SEL_PARSER.parse("resid 1 or not resid 2");
    ASSERT_TRUE(std::dynamic_pointer_cast<OrSelection>(sel_tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<NotSelection>(sel_tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->right->left));
    selected = evaluate_sel_tree(sel_tree, data);
    EXPECT_THAT(*selected, UnorderedElementsAre(0, 1, 2, 3, 6, 7, 8, 9));

    sel_tree = SEL_PARSER.parse("not resid 1 and resid 2");
    ASSERT_TRUE(std::dynamic_pointer_cast<AndSelection>(sel_tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<NotSelection>(sel_tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->left->left));
    selected = evaluate_sel_tree(sel_tree, data);
    EXPECT_THAT(*selected, UnorderedElementsAre(4, 5));

    sel_tree = SEL_PARSER.parse("resid 1 and not resid 2");
    ASSERT_TRUE(std::dynamic_pointer_cast<AndSelection>(sel_tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<NotSelection>(sel_tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->right->left));
    selected = evaluate_sel_tree(sel_tree, data);
    EXPECT_THAT(*selected, UnorderedElementsAre(2, 3));

    sel_tree = SEL_PARSER.parse("resid 1 or not (resid 2 or resid 3)");
    ASSERT_TRUE(std::dynamic_pointer_cast<OrSelection>(sel_tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<NotSelection>(sel_tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<OrSelection>(sel_tree->right->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->right->left->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->right->left->right));
    selected = evaluate_sel_tree(sel_tree, data);
    EXPECT_THAT(*selected, UnorderedElementsAre(0, 1, 2, 3, 8, 9));

    sel_tree = SEL_PARSER.parse("resid 1 and not (resid 2 or resid 3)");
    ASSERT_TRUE(std::dynamic_pointer_cast<AndSelection>(sel_tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<NotSelection>(sel_tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<OrSelection>(sel_tree->right->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->right->left->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->right->left->right));
    selected = evaluate_sel_tree(sel_tree, data);
    EXPECT_THAT(*selected, UnorderedElementsAre(2, 3));
}

TEST(Selection, NumPropParsing) {
    // Resid
    auto sel_tree = SEL_PARSER.parse("resid 0 or resid 2:4 or resid 6:7 9 11 or resid 30");
    // Nodes construction
    ASSERT_TRUE(sel_tree);
    ASSERT_TRUE(sel_tree->left);
    ASSERT_TRUE(sel_tree->right);
    EXPECT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->left->left->left));
    EXPECT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->left->left->right));
    EXPECT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->left->right));
    EXPECT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->right));
    // Selection
    MolData data = create_moldata(12, 1, 1, 1, 1);
    std::shared_ptr<std::unordered_set<index_t>> selected = evaluate_sel_tree(sel_tree, data);
    EXPECT_THAT(*selected, UnorderedElementsAre(0, 2, 3, 4, 6, 7, 9, 11));
}

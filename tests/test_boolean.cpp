#include "selections/boolean.hpp"
#include "selections/SelectionStack.hpp"
#include "selections/SelectionNode.hpp"
#include "selections/SelectionIndices.hpp"
#include <molpp/internal/MolData.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace mol::internal;
using namespace testing;

//! Helper class for testing SelectionNode implementations.
class HelperNode : public SelectionNode
{
public:
    HelperNode(std::unordered_set<index_t> const& selected)
    : selected{selected}
    {}

    void evaluate(SelectionStack& stack, MolData const& /*data*/, Frame /*frame*/) const override
    {
        SelectionIndices indices = stack.pop_indices();
        for (index_t const index : *indices.available)
        {
            if (selected.contains(index))
            {
                indices.selected->insert(index);
            }
        }
    }

    std::unordered_set<index_t> selected;
};

//! Test fixture for boolean selection nodes.
class BooleanTest : public testing::Test
{
public:
    BooleanTest()
    : all{0, 1, 2}
    , left{0, 1}
    , right{1, 2}
    {
        *indices.available = all;
    }

    void evaluate()
    {
        MolData data(0);
        while (!stack.empty_nodes())
        {
            std::shared_ptr<SelectionNode> node = stack.pop_node();
            node->evaluate(stack, data, 0);
        }
    }

    SelectionStack stack;
    std::unordered_set<index_t> all;
    std::unordered_set<index_t> left;
    std::unordered_set<index_t> right;
    SelectionIndices indices;
};

// OrSelection should select the union of the two operands.
TEST_F(BooleanTest, or_selection)
{
    std::shared_ptr<SelectionNode> or_node = std::make_shared<OrSelection>();
    or_node->left = std::make_shared<HelperNode>(left);
    or_node->right = std::make_shared<HelperNode>(right);

    stack.push_indices(indices);
    stack.push_node(or_node);
    evaluate();

    EXPECT_THAT(*indices.selected, UnorderedElementsAre(0, 1, 2));
}

// AndSelection should select the intersection of the two operands.
TEST_F(BooleanTest, and_selection)
{
    std::shared_ptr<SelectionNode> and_node = std::make_shared<AndSelection>();
    and_node->left = std::make_shared<HelperNode>(left);
    and_node->right = std::make_shared<HelperNode>(right);

    stack.push_indices(indices);
    stack.push_node(and_node);
    evaluate();

    EXPECT_THAT(*indices.selected, UnorderedElementsAre(1));
}

// NotSelection should select the inverse of the operand.
TEST_F(BooleanTest, not_selection)
{
    std::shared_ptr<SelectionNode> not_node = std::make_shared<NotSelection>();
    not_node->left = std::make_shared<HelperNode>(left);

    stack.push_indices(indices);
    stack.push_node(not_node);
    evaluate();

    EXPECT_THAT(*indices.selected, UnorderedElementsAre(2));
}

// AllSelection should select all available indices.
TEST_F(BooleanTest, all_selection)
{
    std::shared_ptr<SelectionNode> all_node = std::make_shared<AllSelection>();

    stack.push_indices(indices);
    stack.push_node(all_node);
    evaluate();

    EXPECT_THAT(*indices.selected, UnorderedElementsAre(0, 1, 2));
}

#include "selections/SelectionStack.hpp"
#include "selections/SelectionNode.hpp"
#include "selections/SelectionIndices.hpp"

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace mol::internal;
using namespace testing;

namespace mol::internal
{
class MolData;
class SelectionStack;
} // namespace mol::internal

//! Dummy SelectionNode for testing purposes.
class DummyNode : public SelectionNode
{
public:
    void evaluate(SelectionStack& stack, MolData const& data, Frame frame) const override
    {}
};

//! Test fixture for SelectionStack.
class SelectionStackTest : public testing::Test
{
public:
    SelectionStackTest()
    : node1{std::make_shared<DummyNode>()}
    , node2{std::make_shared<DummyNode>()}
    {}

    SelectionStack stack;
    std::shared_ptr<SelectionNode> node1;
    std::shared_ptr<SelectionNode> node2;
    SelectionIndices indices1;
    SelectionIndices indices2;
};

// SelectionStack should be empty by default.
TEST_F(SelectionStackTest, empty_by_default)
{
    EXPECT_TRUE(stack.empty_nodes());
    EXPECT_TRUE(stack.empty_indices());
}

// SelectionStack nodes should not be empty after pushing a node1.
TEST_F(SelectionStackTest, not_empty_after_pushing_node)
{
    stack.push_node(node1);

    EXPECT_FALSE(stack.empty_nodes());
}

// SelectionStack indices1 should not be empty after pushing indices1.
TEST_F(SelectionStackTest, not_empty_after_pushing_indices)
{
    stack.push_indices(indices1);

    EXPECT_FALSE(stack.empty_indices());
}

// SelectionStack should be empty after pushing and calling clear.
TEST_F(SelectionStackTest, empty_after_clear)
{
    stack.push_node(node1);
    stack.push_indices(indices1);

    stack.clear();

    EXPECT_TRUE(stack.empty_nodes());
    EXPECT_TRUE(stack.empty_indices());
}

// SelectionStack::pop_node should return pushed nodes from last to first.
TEST_F(SelectionStackTest, pop_node)
{
    stack.push_node(node1);
    stack.push_node(node2);

    EXPECT_EQ(stack.pop_node(), node2);
    EXPECT_EQ(stack.pop_node(), node1);
}

// SelectionStack::pop_indices should return pushed indices from last to first.
TEST_F(SelectionStackTest, pop_indices)
{
    stack.push_indices(indices1);
    stack.push_indices(indices2);

    EXPECT_EQ(stack.pop_indices(), indices2);
    EXPECT_EQ(stack.pop_indices(), indices1);
}

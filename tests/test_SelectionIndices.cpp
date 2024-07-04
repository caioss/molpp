#include "selections/SelectionIndices.hpp"

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace mol::internal;
using namespace testing;

//! Test fixture for SelectionIndices.
class SelectionIndicesTest : public testing::Test
{
public:
    SelectionIndicesTest()
    : available{std::make_shared<std::unordered_set<index_t>>()}
    , selected{std::make_shared<std::unordered_set<index_t>>()}
    {}

    std::shared_ptr<std::unordered_set<index_t>> available;
    std::shared_ptr<std::unordered_set<index_t>> selected;
};

// Default constructor should create empty available and selected sets
TEST_F(SelectionIndicesTest, default_constructor)
{
    SelectionIndices indices;

    ASSERT_TRUE(indices.available);
    ASSERT_TRUE(indices.selected);
    EXPECT_TRUE(indices.available->empty());
    EXPECT_TRUE(indices.selected->empty());
}

// Constructor with available and selected sets should set them accordingly
TEST_F(SelectionIndicesTest, constructor_with_sets)
{
    SelectionIndices indices(available, selected);

    EXPECT_EQ(indices.available, available);
    EXPECT_EQ(indices.selected, selected);
}

// Equality operator should return true for equal SelectionIndices
TEST_F(SelectionIndicesTest, equality_operator)
{
    SelectionIndices indices1(available, selected);
    SelectionIndices indices2(available, selected);

    EXPECT_EQ(indices1, indices2);
}

// Equality operator should return false for different available sets
TEST_F(SelectionIndicesTest, equality_operator_available)
{
    std::shared_ptr<std::unordered_set<index_t>> available2 = std::make_shared<std::unordered_set<index_t>>();

    SelectionIndices indices1(available, selected);
    SelectionIndices indices2(available2, selected);

    EXPECT_NE(indices1, indices2);
}

// Equality operator should return false for different selected sets
TEST_F(SelectionIndicesTest, equality_operator_selected)
{
    std::shared_ptr<std::unordered_set<index_t>> available = std::make_shared<std::unordered_set<index_t>>();
    std::shared_ptr<std::unordered_set<index_t>> selected2 = std::make_shared<std::unordered_set<index_t>>();

    SelectionIndices indices1(available, selected);
    SelectionIndices indices2(available, selected2);

    EXPECT_NE(indices1, indices2);
}

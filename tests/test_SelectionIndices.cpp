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
    : available{SelectionIndices::make_index_set()}
    , selected{SelectionIndices::make_index_set()}
    {}

    SelectionIndices::IndexSet available;
    SelectionIndices::IndexSet selected;
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

// make_index_set should return a new empty indices set
TEST_F(SelectionIndicesTest, make_index_set)
{
    SelectionIndices::IndexSet indices1 = SelectionIndices::make_index_set();
    SelectionIndices::IndexSet indices2 = SelectionIndices::make_index_set();

    ASSERT_TRUE(indices1);
    ASSERT_TRUE(indices2);
    EXPECT_TRUE(indices1->empty());
    EXPECT_TRUE(indices2->empty());
    EXPECT_NE(indices1, indices2);
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
    SelectionIndices::IndexSet available2 = SelectionIndices::make_index_set();

    SelectionIndices indices1(available, selected);
    SelectionIndices indices2(available2, selected);

    EXPECT_NE(indices1, indices2);
}

// Equality operator should return false for different selected sets
TEST_F(SelectionIndicesTest, equality_operator_selected)
{
    SelectionIndices::IndexSet available = SelectionIndices::make_index_set();
    SelectionIndices::IndexSet selected2 = SelectionIndices::make_index_set();

    SelectionIndices indices1(available, selected);
    SelectionIndices indices2(available, selected2);

    EXPECT_NE(indices1, indices2);
}

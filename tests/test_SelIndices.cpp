#include "utils.hpp"

#include <molpp/internal/SelIndices.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <array>
#include <ranges>
#include <iterator>

using namespace testing;

class SelIndicesTest : public ::testing::Test
{
public:
    SelIndicesTest()
    : indices_from_list(std::to_array<mol::index_t>({1, 2, 3, 4, 5}), false)
    , indices_as_is(std::to_array<mol::index_t>({2, 2, 0, 1, 1}), true)
    , indices_from_size{5}
    {
    }

    mol::internal::SelIndices indices_from_list;
    mol::internal::SelIndices indices_as_is;
    mol::internal::SelIndices indices_from_size;
};

// SelIndices::size should be equal to the number of indices
TEST_F(SelIndicesTest, size)
{
    EXPECT_EQ(indices_from_list.size(), 5);
    EXPECT_EQ(indices_as_is.size(), 5);
    EXPECT_EQ(indices_from_size.size(), 5);
}

// SelIndices::indices should return all indices
TEST_F(SelIndicesTest, indices)
{
    EXPECT_THAT(indices_from_list.indices(), ElementsAre(1, 2, 3, 4, 5));
    EXPECT_THAT(indices_as_is.indices(), ElementsAre(2, 2, 0, 1, 1));
    EXPECT_THAT(indices_from_size.indices(), ElementsAre(0, 1, 2, 3, 4));
}

// Constructor with keep set to false should remove duplicates and sort them
TEST_F(SelIndicesTest, remove_duplicates_and_sort)
{
    mol::internal::SelIndices indices(std::to_array<mol::index_t>({2, 1, 3, 4, 5, 5, 4, 3, 2, 1}), false);
    EXPECT_THAT(indices.indices(), ElementsAre(1, 2, 3, 4, 5));
}

// Constructor with keep equals to true should keep indices as is
TEST_F(SelIndicesTest, keep_indices)
{
    mol::internal::SelIndices indices(std::to_array<mol::index_t>({2, 1, 3, 4, 5, 5, 4, 3, 2, 1}), true);
    EXPECT_THAT(indices.indices(), ElementsAre(2, 1, 3, 4, 5, 5, 4, 3, 2, 1));
}

// SelIndices::contains should return true if the index is in the list
TEST_F(SelIndicesTest, contains)
{
    EXPECT_TRUE(indices_from_list.contains(1));
    EXPECT_TRUE(indices_from_list.contains(2));
    EXPECT_TRUE(indices_from_list.contains(3));
    EXPECT_TRUE(indices_from_list.contains(4));
    EXPECT_TRUE(indices_from_list.contains(5));

    EXPECT_TRUE(indices_as_is.contains(0));
    EXPECT_TRUE(indices_as_is.contains(1));
    EXPECT_TRUE(indices_as_is.contains(2));

    EXPECT_TRUE(indices_from_size.contains(0));
    EXPECT_TRUE(indices_from_size.contains(1));
    EXPECT_TRUE(indices_from_size.contains(2));
    EXPECT_TRUE(indices_from_size.contains(3));
    EXPECT_TRUE(indices_from_size.contains(4));
}

// SelIndices::contains should return false if the index is not in the list
TEST_F(SelIndicesTest, does_not_contain)
{
    EXPECT_FALSE(indices_from_list.contains(0));
    EXPECT_FALSE(indices_from_list.contains(6));

    EXPECT_FALSE(indices_as_is.contains(3));

    EXPECT_FALSE(indices_from_size.contains(5));
}

// SelIndices::begin should return a const iterator to the beginning of the list
TEST_F(SelIndicesTest, begin)
{
    auto it = indices_from_list.begin();

    EXPECT_EQ(*it, 1);
}

// SelIndices::end should return a const iterator to the end of the list
TEST_F(SelIndicesTest, end)
{
    auto it = indices_from_list.end();

    EXPECT_EQ(*std::prev(it), 5);
}

// SelIndices::begin and SelIndices::end should form a valid range
TEST_F(SelIndicesTest, valid_range)
{
    auto begin = indices_from_list.begin();
    auto end = indices_from_list.end();
    auto range = std::ranges::subrange(begin, end);

    EXPECT_THAT(view2vector(range), ElementsAre(1, 2, 3, 4, 5));
}

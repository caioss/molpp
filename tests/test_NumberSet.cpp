#include "selections/NumberSet.hpp"

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace mol::internal;
using namespace testing;

//! Test fixture for boolean selection nodes.
class NumberSetTest : public testing::Test
{
public:
    NumberSetTest()
    {
        numbers.add_number({"0"});
        numbers.add_number({"2"});
        numbers.add_range({"4", "6"});
        numbers.add_range({"8", "10"});
    }

    NumberSet numbers;
};

// NumberSet::has should return true for numbers in the set
TEST_F(NumberSetTest, has_number)
{
    EXPECT_TRUE(numbers.has(0));
    EXPECT_TRUE(numbers.has(2));
    EXPECT_TRUE(numbers.has(4));
    EXPECT_TRUE(numbers.has(5));
    EXPECT_TRUE(numbers.has(6));
    EXPECT_TRUE(numbers.has(8));
    EXPECT_TRUE(numbers.has(9));
    EXPECT_TRUE(numbers.has(10));
}

// NumberSet::has should return false for numbers not in the set
TEST_F(NumberSetTest, has_not_number)
{
    EXPECT_FALSE(numbers.has(-1));
    EXPECT_FALSE(numbers.has(1));
    EXPECT_FALSE(numbers.has(3));
    EXPECT_FALSE(numbers.has(7));
    EXPECT_FALSE(numbers.has(11));
}

// Adding a number to the set should make NumberSet::has return true for that number
TEST_F(NumberSetTest, add_number)
{
    numbers.add_number({"1"});
    EXPECT_TRUE(numbers.has(1));
}

// Adding a range to the set should make NumberSet::has return true for numbers in that range
TEST_F(NumberSetTest, add_range)
{
    numbers.add_range({"12", "14"});
    EXPECT_TRUE(numbers.has(12));
    EXPECT_TRUE(numbers.has(13));
    EXPECT_TRUE(numbers.has(14));
}

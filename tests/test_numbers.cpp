#include "selections/numbers.hpp"

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace mol::internal;
using namespace testing;

// Test floating point number comparison.
TEST(SelNumberTest, floating_point)
{
    SelNumber number("1.0001");

    EXPECT_EQ(number, 1.0001);
    EXPECT_EQ(number, 1.0002);
    EXPECT_NE(number, 1.0003);
}

// Test integer number comparison.
TEST(SelNumberTest, integer)
{
    SelNumber number("1");

    EXPECT_EQ(number, 1);
    EXPECT_NE(number, 2);
}

// Test range comparison to integers.
TEST(SelNumberRangeTest, integer)
{
    SelNumberRange range("3", "5");

    EXPECT_FALSE(range.has(-1));
    EXPECT_FALSE(range.has(2));
    EXPECT_TRUE(range.has(3));
    EXPECT_TRUE(range.has(4));
    EXPECT_TRUE(range.has(5));
    EXPECT_FALSE(range.has(6));
}

// Test range comparison to floating point numbers.
TEST(SelNumberRangeTest, floating_point)
{
    SelNumberRange range("3", "5");

    EXPECT_FALSE(range.has(-1.0));
    EXPECT_FALSE(range.has(2.0));
    EXPECT_TRUE(range.has(3.0));
    EXPECT_TRUE(range.has(4.0));
    EXPECT_TRUE(range.has(5.0));
    EXPECT_FALSE(range.has(6.0));
}

// Test range with consecutive limits.
TEST(SelNumberRangeTest, consecutive)
{
    SelNumberRange range("3", "4");

    EXPECT_FALSE(range.has(2));
    EXPECT_TRUE(range.has(3));
    EXPECT_TRUE(range.has(4));
    EXPECT_FALSE(range.has(5));
}

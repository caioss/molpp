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

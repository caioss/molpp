#include "tools/math.hpp"

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace testing;

// approximately_equal
TEST(Math, approximately_equal)
{
    EXPECT_TRUE(mol::internal::approximately_equal(95.1, 100.0, 0.05));
    EXPECT_FALSE(mol::internal::approximately_equal(95.1, 100.0, 0.01));
}

// essentially_equal
TEST(Math, essentially_equal)
{
    EXPECT_FALSE(mol::internal::essentially_equal(95.1, 100.0, 0.05));
}

// definitely_greater
TEST(Math, definitely_greater)
{
    EXPECT_TRUE(mol::internal::definitely_greater(106.0, 100.0, 0.05));
    EXPECT_FALSE(mol::internal::definitely_greater(95.1, 100.0, 0.05));
}

// definitely_less
TEST(Math, definitely_less)
{
    EXPECT_FALSE(mol::internal::definitely_less(105.1, 100.0, 0.05));
    EXPECT_TRUE(mol::internal::definitely_less(94.9, 100.0, 0.05));
}

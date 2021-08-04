#include "tools/iterators.hpp"
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <vector>

using namespace mol::internal;
using namespace testing;

TEST(Iterators, IteratorWrapper) {
    using iter_type = typename std::vector<int>::iterator;
    std::vector<int> values(2);

    class Wrapper : public IteratorWrapper<iter_type>
    {
        using base = IteratorWrapper<iter_type>;
        using base::IteratorWrapper;
    };

    Wrapper iter1(values.begin());
    Wrapper iter2(values.begin());
    Wrapper end(values.end());

    EXPECT_EQ(iter1, iter2);
    EXPECT_EQ(iter1++, iter2);
    EXPECT_TRUE(iter1 != iter2);
    EXPECT_EQ(++iter1, end);
    EXPECT_EQ(iter1 - iter2, 2);
}

TEST(Iterators, Range) {
    using iter_type = typename std::vector<int>::iterator;
    using RangeInt = mol::internal::Range<iter_type>;
    std::vector<int> values(1);

    RangeInt range(values.begin(), values.end());
    EXPECT_EQ(range.begin(), values.begin());
    EXPECT_EQ(range.end(), values.end());
    EXPECT_TRUE(range.is_valid());

    RangeInt invalid(values.end(), values.end());
    EXPECT_EQ(invalid.begin(), values.end());
    EXPECT_EQ(invalid.end(), values.end());
    EXPECT_FALSE(invalid.is_valid());
}

#include "utils.hpp"
#include "tools/algorithms.hpp"
#include "tools/math.hpp"
#include <molpp/Common.hpp>
#include <molpp/internal/SelIndices.hpp>
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <vector>
#include <numeric>

using namespace mol::internal;
using namespace testing;

TEST(Math, Comparison) {
    EXPECT_TRUE(approximately_equal(95.1, 100.0, 0.05));
    EXPECT_FALSE(essentially_equal(95.1, 100.0, 0.05));

    EXPECT_TRUE(definitely_greater(106.0, 100.0, 0.05));
    EXPECT_FALSE(definitely_greater(95.1, 100.0, 0.05));

    EXPECT_FALSE(definitely_less(105.1, 100.0, 0.05));
    EXPECT_TRUE(definitely_less(94.9, 100.0, 0.05));
}

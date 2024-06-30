#include "auxiliary.hpp"

#include <molpp/internal/SegmentData.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol::internal;
using namespace testing;

//! Test fixture for SegmentData class
class SegmentDataTest : public ::testing::Test
{
public:
    SegmentDataTest()
    {
        data.resize(1);
    }

    SegmentData data;
};

// Size
TEST_F(SegmentDataTest, size)
{
    EXPECT_EQ(data.size(), 1);
}

// Resize should change size
TEST_F(SegmentDataTest, resize)
{
    data.resize(2);

    EXPECT_EQ(data.size(), 2);
}

// Default name should be an empty string
TEST_F(SegmentDataTest, default_name)
{
    EXPECT_EQ(data.name(0), "");
}

// Name should be settable
TEST_F(SegmentDataTest, set_name)
{
    ASSERT_EQ(data.name(0), "");

    data.name(0) = "A";

    EXPECT_EQ(data.name(0), "A");
}

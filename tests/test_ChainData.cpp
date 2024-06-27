#include "auxiliary.hpp"

#include <molpp/internal/ChainData.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol::internal;
using namespace testing;

//! Test fixture for ChainData class
class ChainDataTest : public ::testing::Test
{
public:
    ChainDataTest()
    {
        data.resize(1);
    }

    ChainData data;
};

// Size
TEST_F(ChainDataTest, size)
{
    EXPECT_EQ(data.size(), 1);
}

// Resize should change size
TEST_F(ChainDataTest, resize)
{
    data.resize(2);

    EXPECT_EQ(data.size(), 2);
}

// Default name should be an empty string
TEST_F(ChainDataTest, default_name)
{
    EXPECT_EQ(data.name(0), "");
}

// Name should be settable
TEST_F(ChainDataTest, set_name)
{
    ASSERT_EQ(data.name(0), "");

    data.name(0) = "A";

    EXPECT_EQ(data.name(0), "A");
}

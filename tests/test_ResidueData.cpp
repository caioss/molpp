#include "auxiliary.hpp"

#include <molpp/internal/ResidueData.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol::internal;
using namespace testing;

//! Test fixture for ResidueData class
class ResidueDataTest : public ::testing::Test
{
public:
    ResidueDataTest()
    {
        data.resize(1);
    }

    ResidueData data;
};

// Size
TEST_F(ResidueDataTest, size)
{
    EXPECT_EQ(data.size(), 1);
}

// Resize should change size
TEST_F(ResidueDataTest, resize)
{
    data.resize(2);

    EXPECT_EQ(data.size(), 2);
}

// Default ID should be -1
TEST_F(ResidueDataTest, default_id)
{
    EXPECT_EQ(data.id(0), -1);
}

// Default name should be an empty string
TEST_F(ResidueDataTest, default_name)
{
    EXPECT_EQ(data.name(0), "");
}

// Default segid should be an empty string
TEST_F(ResidueDataTest, default_segid)
{
    EXPECT_EQ(data.segid(0), "");
}

// Default chain should be an empty string
TEST_F(ResidueDataTest, default_chain)
{
    EXPECT_EQ(data.chain(0), "");
}

// ResidueData::set should set all data of a residue
TEST_F(ResidueDataTest, set)
{
    data.set(0, 1, "A", "B", "C");

    EXPECT_EQ(data.id(0), 1);
    EXPECT_EQ(data.name(0), "A");
    EXPECT_EQ(data.segid(0), "B");
    EXPECT_EQ(data.chain(0), "C");
}

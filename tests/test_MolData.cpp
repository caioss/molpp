#include "utils.hpp"

#include <molpp/internal/MolData.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace testing;

//! Test fixture for MolData
class MolDataTest : public ::testing::Test
{
public:
    MolDataTest()
    : data(create_moldata(4, 1, 4, 4, 2))
    {
    }

    mol::internal::MolData data;
};

// MolData::size should return the correct number of atoms
TEST_F(MolDataTest, size_atoms)
{
    EXPECT_EQ(data.size<mol::Atom>(), 4);
}

// MolData::size should return the correct number of residues
TEST_F(MolDataTest, size_residues)
{
    EXPECT_EQ(data.size<mol::Residue>(), 4);
}

// MolData::size should return the correct number of chains
TEST_F(MolDataTest, size_chains)
{
    EXPECT_EQ(data.size<mol::Chain>(), 4);
}

// MolData::size should return the correct number of segments
TEST_F(MolDataTest, size_segments)
{
    EXPECT_EQ(data.size<mol::Segment>(), 4);
}

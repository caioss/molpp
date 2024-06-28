#include "auxiliary.hpp"

#include <molpp/internal/MolData.hpp>
#include <molpp/Chain.hpp>
#include <molpp/Residue.hpp>
#include <molpp/ResidueSel.hpp>
#include <molpp/MolError.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace testing;

//! Test fixture for Chain class
class ChainTest : public ::testing::Test
{
public:
    ChainTest()
    : data(create_moldata(3, 1, 3, 1, 1))
    , chain(1, 0, data)
    , const_chain(1, 0, data)
    , null_frame_chain(0, std::nullopt, data)
    {}

    mol::internal::MolData data;
    Chain chain;
    Chain const const_chain;
    Chain null_frame_chain;
};

// Chain::category should return the correct category
TEST_F(ChainTest, category)
{
    EXPECT_EQ(Chain::category(), MolecularEntityCategory::Chain);
}

// Chain::size should return the number of residues in the chain
TEST_F(ChainTest, size)
{
    EXPECT_EQ(chain.size(), 1);
}

// Chain::index should return the correct index
TEST_F(ChainTest, index)
{
    EXPECT_EQ(chain.index(), 1);
    EXPECT_EQ(const_chain.index(), 1);
}

// Chain::indices should return a list with only the index
TEST_F(ChainTest, indices)
{
    EXPECT_THAT(chain.indices(), ElementsAre(1));
    EXPECT_THAT(const_chain.indices(), ElementsAre(1));
}

// Chain::frame should return the correct frame
TEST_F(ChainTest, frames)
{
    EXPECT_EQ(chain.frame(), 0);
    EXPECT_EQ(const_chain.frame(), 0);
    EXPECT_FALSE(null_frame_chain.frame());
}

// Setting a valid frame should update the frame
TEST_F(ChainTest, set_valid_frame)
{
    chain.set_frame(0);
    EXPECT_EQ(chain.frame(), 0);
}

// Setting a frame to nullopt should update the frame
TEST_F(ChainTest, set_null_frame)
{
    chain.set_frame(std::nullopt);
    EXPECT_FALSE(chain.frame());
}

// Setting an invalid frame should throw a MolError
TEST_F(ChainTest, set_invalid_frame)
{
    EXPECT_THROW(chain.set_frame(1), MolError);
}

// Chain::name should return the correct name
TEST_F(ChainTest, name)
{
    EXPECT_EQ(chain.name(), "B");
    EXPECT_EQ(const_chain.name(), "B");
}

// Chain::set_name should update the name
TEST_F(ChainTest, set_name)
{
    chain.set_name("A");

    EXPECT_EQ(chain.name(), "A");
}

// Chain::add_residue with an index should make the residue exclusive to the new chain
TEST_F(ChainTest, add_residue_from_index)
{
    Residue new_residue(0, std::nullopt, data);
    Chain old_chain = new_residue.chain().value();
    ASSERT_NE(old_chain, chain);

    chain.add_residue(new_residue.index());
    ResidueSel chain_residues(chain);
    ResidueSel old_chain_residues(old_chain);

    EXPECT_THAT(chain_residues.indices(), UnorderedElementsAre(0, 1));
    EXPECT_THAT(old_chain_residues.indices(), ElementsAre());
    EXPECT_EQ(new_residue.chain_index(), chain.index());
    EXPECT_EQ(chain.size(), 2);
    EXPECT_EQ(old_chain.size(), 0);
}

// Chain::add_residue with an Residue should make the residue exclusive to the new chain
TEST_F(ChainTest, add_residue_from_residue)
{
    Residue new_residue(0, std::nullopt, data);
    Chain old_chain = new_residue.chain().value();
    ASSERT_NE(old_chain, chain);

    chain.add_residue(new_residue);
    ResidueSel chain_residues(chain);
    ResidueSel old_chain_residues(old_chain);

    EXPECT_THAT(chain_residues.indices(), UnorderedElementsAre(0, 1));
    EXPECT_THAT(old_chain_residues.indices(), ElementsAre());
    EXPECT_EQ(new_residue.chain_index(), chain.index());
    EXPECT_EQ(chain.size(), 2);
    EXPECT_EQ(old_chain.size(), 0);
}

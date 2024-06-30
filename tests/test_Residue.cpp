#include "auxiliary.hpp"

#include <molpp/internal/MolData.hpp>
#include <molpp/Residue.hpp>
#include <molpp/Atom.hpp>
#include <molpp/Chain.hpp>
#include <molpp/Segment.hpp>
#include <molpp/AtomSel.hpp>
#include <molpp/MolError.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace testing;

//! Test fixture for Residue class
class ResidueTest : public ::testing::Test
{
public:
    ResidueTest()
    : data(create_moldata(3, 1, 2, 2, 1))
    , residue(1, 0, data)
    , const_residue(1, 0, data)
    , null_frame_residue(0, std::nullopt, data)
    {}

    mol::internal::MolData data;
    Residue residue;
    Residue const const_residue;
    Residue null_frame_residue;
};

// Residue::category should return the correct category
TEST_F(ResidueTest, category)
{
    EXPECT_EQ(Residue::category(), MolecularEntityCategory::Residue);
}

// Residue::index should return the correct index
TEST_F(ResidueTest, index)
{
    EXPECT_EQ(residue.index(), 1);
    EXPECT_EQ(const_residue.index(), 1);
}

// Residue::indices should return a list with only the index
TEST_F(ResidueTest, indices)
{
    EXPECT_THAT(residue.indices(), ElementsAre(1));
    EXPECT_THAT(const_residue.indices(), ElementsAre(1));
}

// Residue::frames should return the correct frame
TEST_F(ResidueTest, frames)
{
    EXPECT_EQ(residue.frame(), 0);
    EXPECT_EQ(const_residue.frame(), 0);
    EXPECT_FALSE(null_frame_residue.frame());
}

// Setting a valid frame should update the frame
TEST_F(ResidueTest, set_valid_frame)
{
    residue.set_frame(0);
    EXPECT_EQ(residue.frame(), 0);
}

// Setting a null frame should update the frame
TEST_F(ResidueTest, set_null_frame)
{
    residue.set_frame(std::nullopt);
    EXPECT_FALSE(residue.frame());
}

// Setting an invalid frame should throw a MolError
TEST_F(ResidueTest, set_invalid_frame)
{
    EXPECT_THROW(residue.set_frame(1), MolError);
}

// Residue::chain should return the correct chain
TEST_F(ResidueTest, chain)
{
    ASSERT_TRUE(residue.chain());
    EXPECT_EQ(residue.chain(), Chain(1, residue.frame(), data));
}

// Residue::chain_index should return the correct chain index
TEST_F(ResidueTest, chain_index)
{
    ASSERT_TRUE(residue.chain_index());
    EXPECT_EQ(residue.chain_index(), 1);
    EXPECT_EQ(const_residue.chain_index(), 1);
}

// Residue::segment should return the correct segment
TEST_F(ResidueTest, segment)
{
    ASSERT_TRUE(residue.segment());
    EXPECT_EQ(residue.segment(), Segment(1, residue.frame(), data));
}

// Residue::segment_index should return the correct segment index
TEST_F(ResidueTest, segment_index)
{
    ASSERT_TRUE(residue.segment_index());
    EXPECT_EQ(residue.segment_index(), 1);
    EXPECT_EQ(const_residue.segment_index(), 1);
}

// Residue::name should return the correct name
TEST_F(ResidueTest, name_property)
{
    EXPECT_EQ(const_residue.name(), "B");
    EXPECT_EQ(residue.name(), "B");
}

// Residue::set_name should update the name
TEST_F(ResidueTest, set_name_property)
{
    residue.set_name("C");

    EXPECT_EQ(residue.name(), "C");
}

// Residue::id should return the correct ID
TEST_F(ResidueTest, id_property)
{
    EXPECT_EQ(const_residue.id(), 1);
    EXPECT_EQ(residue.id(), 1);
}

// Residue::set_id should update the ID
TEST_F(ResidueTest, set_id_property)
{
    residue.set_id(2);

    EXPECT_EQ(residue.id(), 2);
}

// Residue::size should return the number of atoms in the residue
TEST_F(ResidueTest, size)
{
    EXPECT_EQ(residue.size(), 1);
}

// Residue::add_atom with an index should make the atom exclusive to the new residue
TEST_F(ResidueTest, add_atom_from_index)
{
    Atom new_atom(0, std::nullopt, data);
    Residue old_residue = new_atom.residue().value();
    ASSERT_NE(old_residue, residue);

    residue.add_atom(new_atom.index());
    AtomSel residue_atoms(residue);
    AtomSel old_residue_atoms(old_residue);

    EXPECT_THAT(residue_atoms.indices(), UnorderedElementsAre(0, 1));
    EXPECT_THAT(old_residue_atoms.indices(), ElementsAre());
    EXPECT_EQ(new_atom.residue_index(), residue.index());
    EXPECT_EQ(residue.size(), 2);
    EXPECT_EQ(old_residue.size(), 0);
}

// Residue::add_atom with an Atom should make the atom exclusive to the new residue
TEST_F(ResidueTest, add_atom_from_atom)
{
    Atom new_atom(0, std::nullopt, data);
    Residue old_residue = new_atom.residue().value();
    ASSERT_NE(old_residue, residue);

    residue.add_atom(new_atom);
    AtomSel residue_atoms(residue);
    AtomSel old_residue_atoms(old_residue);

    EXPECT_THAT(residue_atoms.indices(), UnorderedElementsAre(0, 1));
    EXPECT_THAT(old_residue_atoms.indices(), ElementsAre());
    EXPECT_EQ(new_atom.residue_index(), residue.index());
    EXPECT_EQ(residue.size(), 2);
    EXPECT_EQ(old_residue.size(), 0);
}

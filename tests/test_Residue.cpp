#include "auxiliary.hpp"

#include <molpp/internal/MolData.hpp>
#include <molpp/Residue.hpp>
#include <molpp/Atom.hpp>
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
    : data(create_moldata(3, 1, 1, 1, 1))
    , residue(1, 0, data)
    , const_residue(1, 0, data)
    , null_frame_residue(0, std::nullopt, data)
    {}

    mol::internal::MolData data;
    Residue residue;
    Residue const const_residue;
    Residue null_frame_residue;
};

TEST_F(ResidueTest, category)
{
    EXPECT_EQ(Residue::category(), MolecularEntityCategory::Residue);
}

TEST_F(ResidueTest, index)
{
    EXPECT_EQ(residue.index(), 1);
    EXPECT_EQ(const_residue.index(), 1);
}

TEST_F(ResidueTest, indices)
{
    EXPECT_THAT(residue.indices(), ElementsAre(1));
    EXPECT_THAT(const_residue.indices(), ElementsAre(1));
}

TEST_F(ResidueTest, frames)
{
    EXPECT_EQ(residue.frame(), 0);
    EXPECT_EQ(const_residue.frame(), 0);
    EXPECT_FALSE(null_frame_residue.frame());
}

TEST_F(ResidueTest, set_valid_frame)
{
    residue.set_frame(0);
    EXPECT_EQ(residue.frame(), 0);
}

TEST_F(ResidueTest, set_null_frame)
{
    residue.set_frame(std::nullopt);
    EXPECT_FALSE(residue.frame());
}

TEST_F(ResidueTest, set_invalid_frame)
{
    EXPECT_THROW(residue.set_frame(1), MolError);
}

TEST_F(ResidueTest, residue_id_property)
{
    EXPECT_EQ(const_residue.id(), 1);
    EXPECT_EQ(residue.id(), 1);
}

TEST_F(ResidueTest, set_residue_id_property)
{
    residue.set_id(2);

    EXPECT_EQ(residue.id(), 2);
}

TEST_F(ResidueTest, size)
{
    EXPECT_EQ(residue.size(), 1);
}

// Residue::add_atom with an index should make the atom exclusive to the new residue
TEST_F(ResidueTest, add_atom_from_index)
{
    Residue old_residue(0, std::nullopt, data);
    Atom new_atom(0, std::nullopt, data);

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
    Residue old_residue(0, std::nullopt, data);
    Atom new_atom(0, std::nullopt, data);

    residue.add_atom(new_atom);
    AtomSel residue_atoms(residue);
    AtomSel old_residue_atoms(old_residue);

    EXPECT_THAT(residue_atoms.indices(), UnorderedElementsAre(0, 1));
    EXPECT_THAT(old_residue_atoms.indices(), ElementsAre());
    EXPECT_EQ(new_atom.residue_index(), residue.index());
    EXPECT_EQ(residue.size(), 2);
    EXPECT_EQ(old_residue.size(), 0);
}

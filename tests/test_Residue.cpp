#include "auxiliary.hpp"

#include <molpp/internal/MolData.hpp>
#include <molpp/Residue.hpp>
#include <molpp/Atom.hpp>
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

    MolData data;
    Residue residue;
    Residue const const_residue;
    Residue null_frame_residue;
};

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
    EXPECT_EQ(const_residue.residue_id(), 1);
    EXPECT_EQ(residue.residue_id(), 1);
}

TEST_F(ResidueTest, set_residue_id_property)
{
    residue.set_residue_id(2);

    EXPECT_EQ(residue.residue_id(), 2);
}

TEST_F(ResidueTest, size)
{
    EXPECT_EQ(residue.size(), 1);
}

TEST_F(ResidueTest, add_atom_from_index)
{
    Residue old_residue(0, std::nullopt, data);
    Atom new_atom(0, std::nullopt, data);

    residue.add_atom(new_atom.index());

    EXPECT_THAT(view2vector(residue.as_atom_indices()), UnorderedElementsAre(0, 1));
    EXPECT_THAT(view2vector(old_residue.as_atom_indices()), ElementsAre());
    EXPECT_EQ(new_atom.residue_id(), residue.index());
    EXPECT_EQ(residue.size(), 2);
    EXPECT_EQ(old_residue.size(), 0);
}

TEST_F(ResidueTest, add_atom_from_atom)
{
    Residue old_residue(0, std::nullopt, data);
    Atom new_atom(0, std::nullopt, data);

    residue.add_atom(new_atom);

    EXPECT_THAT(view2vector(residue.as_atom_indices()), UnorderedElementsAre(0, 1));
    EXPECT_THAT(view2vector(old_residue.as_atom_indices()), ElementsAre());
    EXPECT_EQ(new_atom.residue_id(), residue.index());
    EXPECT_EQ(residue.size(), 2);
    EXPECT_EQ(old_residue.size(), 0);
}

TEST_F(ResidueTest, as_atom_indices)
{
    EXPECT_THAT(view2vector(residue.as_atom_indices()), UnorderedElementsAre(1));
}

#include "auxiliary.hpp"

#include "readers/ResidueDetect.hpp"
#include <molpp/internal/MolData.hpp>
#include <molpp/Atom.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol::internal;
using namespace testing;

//! Test fixture for ResidueDetect class
class ResidueDetectTest : public ::testing::Test
{
public:
    ResidueDetectTest()
    : data{4}
    {}

    MolData data;
    ResidueDetect detect;
};

// Register atom with same residue fields should return same residue index
TEST_F(ResidueDetectTest, register_same_residue)
{
    size_t const index = detect.register_atom(1, "ALA", "AA", "A");

    EXPECT_EQ(detect.register_atom(1, "ALA", "AA", "A"), index);
}

// Register atom with different resid should return new residue index
TEST_F(ResidueDetectTest, register_different_resid)
{
    size_t const index = detect.register_atom(1, "ALA", "AA", "A");

    EXPECT_NE(detect.register_atom(2, "ALA", "AA", "A"), index);
}

// Register atom with different chain should return new residue index
TEST_F(ResidueDetectTest, register_different_chain)
{
    size_t const index = detect.register_atom(1, "ALA", "AA", "A");

    EXPECT_NE(detect.register_atom(1, "ALA", "AA", "B"), index);
}

// Register atom with different resname should return new residue index
TEST_F(ResidueDetectTest, register_different_resname)
{
    size_t const index = detect.register_atom(1, "ALA", "AA", "A");

    EXPECT_NE(detect.register_atom(1, "LYS", "AA", "A"), index);
}

// Register atom with different segid should return new residue index
TEST_F(ResidueDetectTest, register_different_segid)
{
    size_t const index = detect.register_atom(1, "ALA", "AA", "A");

    EXPECT_NE(detect.register_atom(1, "ALA", "BB", "A"), index);
}

// Register atom with same residue fields out of order should return same residue index
TEST_F(ResidueDetectTest, register_residues_out_of_order)
{
    size_t const index1 = detect.register_atom(1, "ALA", "AA", "A");
    size_t const index2 = detect.register_atom(1, "LYS", "AA", "A");
    size_t const index3 = detect.register_atom(1, "LYS", "BB", "A");

    EXPECT_EQ(detect.register_atom(1, "LYS", "AA", "A"), index2);
    EXPECT_EQ(detect.register_atom(1, "LYS", "BB", "A"), index3);
    EXPECT_EQ(detect.register_atom(1, "ALA", "AA", "A"), index1);
}

// Update MolData should set residues data
TEST_F(ResidueDetectTest, update_residue_data)
{
    // Register atoms
    data.atoms().residue(0) = detect.register_atom(1, "ALA", "BB", "A");
    data.atoms().residue(1) = detect.register_atom(2, "LYS", "AA", "C");
    data.atoms().residue(2) = detect.register_atom(2, "LYS", "AA", "C");
    data.atoms().residue(3) = detect.register_atom(1, "ALA", "BB", "A");

    // Apply residues data
    detect.update_residue_data(data);
    ResidueData const& residues_data = data.residues();

    // Check residues data
    EXPECT_EQ(residues_data.size(), 2);
    EXPECT_EQ(residues_data.residue_id(0), 1);
    EXPECT_EQ(residues_data.residue_id(1), 2);
    EXPECT_EQ(residues_data.residue_name(0), "ALA");
    EXPECT_EQ(residues_data.residue_name(1), "LYS");
    EXPECT_EQ(residues_data.chain(0), "A");
    EXPECT_EQ(residues_data.chain(1), "C");
    EXPECT_EQ(residues_data.segid(0), "BB");
    EXPECT_EQ(residues_data.segid(1), "AA");
    EXPECT_THAT(view2vector(residues_data.atom_indices(0)), UnorderedElementsAre(0, 3));
    EXPECT_THAT(view2vector(residues_data.atom_indices(1)), UnorderedElementsAre(1, 2));
}

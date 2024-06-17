#include "auxiliary.hpp"

#include "readers/ResidueDetect.hpp"
#include <molpp/MolppCore.hpp>
#include <molpp/internal/MolData.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace mol::internal;
using namespace testing;

//! Test fixture for ResidueDetect class
class ResidueDetectTest : public ::testing::Test
{
public:
    ResidueDetectTest()
    : data{4}
    , topology{data.topology()}
    , detect{data}
    {}

    MolData data;
    Topology const& topology;
    ResidueDetect detect;
};

// Register atom with same residue fields should render same residue
TEST_F(ResidueDetectTest, register_same_residue)
{
    detect.register_atom(0, 1, "ALA", "AA", "A");
    detect.register_atom(1, 1, "ALA", "AA", "A");
    std::vector<index_t> const residue_0 = topology.all_links({MolecularEntityCategory::Residue, 0}, MolecularEntityCategory::Atom);

    EXPECT_THAT(residue_0, UnorderedElementsAre(0, 1));
}

// Register atom with different resid should render different residue
TEST_F(ResidueDetectTest, register_different_resid)
{
    detect.register_atom(0, 1, "ALA", "AA", "A");
    detect.register_atom(1, 2, "ALA", "AA", "A");
    std::vector<index_t> const residue_0 = topology.all_links({MolecularEntityCategory::Residue, 0}, MolecularEntityCategory::Atom);
    std::vector<index_t> const residue_1 = topology.all_links({MolecularEntityCategory::Residue, 1}, MolecularEntityCategory::Atom);

    EXPECT_THAT(residue_0, UnorderedElementsAre(0));
    EXPECT_THAT(residue_1, UnorderedElementsAre(1));
}

// Register atom with different chain should render different residue
TEST_F(ResidueDetectTest, register_different_chain)
{
    detect.register_atom(0, 1, "ALA", "AA", "A");
    detect.register_atom(1, 1, "ALA", "AA", "B");
    std::vector<index_t> const residue_0 = topology.all_links({MolecularEntityCategory::Residue, 0}, MolecularEntityCategory::Atom);
    std::vector<index_t> const residue_1 = topology.all_links({MolecularEntityCategory::Residue, 1}, MolecularEntityCategory::Atom);

    EXPECT_THAT(residue_0, UnorderedElementsAre(0));
    EXPECT_THAT(residue_1, UnorderedElementsAre(1));
}

// Register atom with different resname should render different residue
TEST_F(ResidueDetectTest, register_different_resname)
{
    detect.register_atom(0, 1, "ALA", "AA", "A");
    detect.register_atom(1, 1, "LYS", "AA", "A");
    std::vector<index_t> const residue_0 = topology.all_links({MolecularEntityCategory::Residue, 0}, MolecularEntityCategory::Atom);
    std::vector<index_t> const residue_1 = topology.all_links({MolecularEntityCategory::Residue, 1}, MolecularEntityCategory::Atom);

    EXPECT_THAT(residue_0, UnorderedElementsAre(0));
    EXPECT_THAT(residue_1, UnorderedElementsAre(1));
}

// Register atom with different segid should render different residue
TEST_F(ResidueDetectTest, register_different_segid)
{
    detect.register_atom(0, 1, "ALA", "AA", "A");
    detect.register_atom(1, 1, "ALA", "BB", "A");
    std::vector<index_t> const residue_0 = topology.all_links({MolecularEntityCategory::Residue, 0}, MolecularEntityCategory::Atom);
    std::vector<index_t> const residue_1 = topology.all_links({MolecularEntityCategory::Residue, 1}, MolecularEntityCategory::Atom);

    EXPECT_THAT(residue_0, UnorderedElementsAre(0));
    EXPECT_THAT(residue_1, UnorderedElementsAre(1));
}

// Register atom with same residue fields out of order should create same residue
TEST_F(ResidueDetectTest, register_residues_out_of_order)
{
    detect.register_atom(0, 1, "ALA", "AA", "A");
    detect.register_atom(1, 1, "LYS", "AA", "A");
    detect.register_atom(2, 1, "LYS", "BB", "A");
    detect.register_atom(3, 1, "LYS", "AA", "A");
    detect.register_atom(4, 1, "LYS", "BB", "A");
    detect.register_atom(5, 1, "ALA", "AA", "A");

    std::vector<index_t> const residue_0 = topology.all_links({MolecularEntityCategory::Residue, 0}, MolecularEntityCategory::Atom);
    std::vector<index_t> const residue_1 = topology.all_links({MolecularEntityCategory::Residue, 1}, MolecularEntityCategory::Atom);
    std::vector<index_t> const residue_2 = topology.all_links({MolecularEntityCategory::Residue, 2}, MolecularEntityCategory::Atom);

    EXPECT_THAT(residue_0, UnorderedElementsAre(0, 5));
    EXPECT_THAT(residue_1, UnorderedElementsAre(1, 3));
    EXPECT_THAT(residue_2, UnorderedElementsAre(2, 4));
}

// Update MolData should set structures data
TEST_F(ResidueDetectTest, update_residue_data)
{
    // Register atoms
    detect.register_atom(0, 1, "ALA", "BB", "A");
    detect.register_atom(1, 2, "LYS", "AA", "C");
    detect.register_atom(2, 2, "LYS", "AA", "C");
    detect.register_atom(3, 1, "ALA", "BB", "A");

    // Update data
    detect.update_residue_data(data);

    // Check residues data
    ResidueData const& residues_data = data.residues();
    EXPECT_EQ(residues_data.size(), 2);
    EXPECT_EQ(residues_data.id(0), 1);
    EXPECT_EQ(residues_data.id(1), 2);
    EXPECT_EQ(residues_data.name(0), "ALA");
    EXPECT_EQ(residues_data.name(1), "LYS");
    EXPECT_EQ(residues_data.chain(0), "A");
    EXPECT_EQ(residues_data.chain(1), "C");
    EXPECT_EQ(residues_data.segid(0), "BB");
    EXPECT_EQ(residues_data.segid(1), "AA");
}

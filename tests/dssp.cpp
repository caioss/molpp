#include "analysis/dssp/DSSP.hpp"
#include "matchers.hpp"
#include "auxiliary.hpp"
#include <molpp/MolError.hpp>
#include <molpp/AtomSel.hpp>
#include <molpp/Residue.hpp>
#include <molpp/Property.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <array>

using namespace mol;
using namespace mol::internal;
using namespace testing;

class SSResidueTest : public ::testing::Test
{
public:
    SSResidueTest()
    : mol_data(create_moldata(3, 4, 3, 2, 2))
    , protein_residue(0, 0, &mol_data)
    , proline_residue(1, 0, &mol_data)
    , non_protein_residue(2, 0, &mol_data)
    {
        Name* names = mol_data.properties().add<Atom, Name>(false);
        names->value(0) = "N";
        names->value(1) = "CA";
        names->value(2) = "C";
        names->value(3) = "O";
        protein_residue.set_resname("ALA");

        // Proline
        names->value(4) = "N";
        names->value(5) = "CA";
        names->value(6) = "C";
        names->value(7) = "O";
        proline_residue.set_resname("PRO");

        // Non-proteic
        names->value(8) = "CA";
        names->value(9) = "N";
        non_protein_residue.set_resname("PRO");
    }

    MolData mol_data;
    Residue protein_residue;
    Residue proline_residue;
    Residue non_protein_residue;
};

TEST_F(SSResidueTest, DefaultConstructor)
{
    SSResidue ss_residue;

    EXPECT_FALSE(ss_residue.is_amino_acid());
    EXPECT_EQ(ss_residue.secondary_structure(), mol::Unknown);
    EXPECT_FALSE(ss_residue.is_proline());
    EXPECT_FALSE(ss_residue.is_chain_break());
    EXPECT_EQ(ss_residue.chain(), "");
    EXPECT_FALSE(ss_residue.N());
    EXPECT_FALSE(ss_residue.CA());
    EXPECT_FALSE(ss_residue.C());
    EXPECT_FALSE(ss_residue.O());
}

TEST_F(SSResidueTest, FromProteinResidue)
{
    SSResidue ss_residue(protein_residue);

    EXPECT_TRUE(ss_residue.is_amino_acid());
    EXPECT_EQ(ss_residue.secondary_structure(), mol::Unknown);
    EXPECT_FALSE(ss_residue.is_proline());
    EXPECT_FALSE(ss_residue.is_chain_break());
    EXPECT_EQ(ss_residue.chain(), "A");
    EXPECT_TRUE(ss_residue.N());
    EXPECT_TRUE(ss_residue.CA());
    EXPECT_TRUE(ss_residue.C());
    EXPECT_TRUE(ss_residue.O());
}

TEST_F(SSResidueTest, FromProlineResidue)
{
    SSResidue ss_residue(proline_residue);

    EXPECT_TRUE(ss_residue.is_amino_acid());
    EXPECT_EQ(ss_residue.secondary_structure(), mol::Unknown);
    EXPECT_TRUE(ss_residue.is_proline());
    EXPECT_FALSE(ss_residue.is_chain_break());
    EXPECT_EQ(ss_residue.chain(), "B");
    EXPECT_TRUE(ss_residue.N());
    EXPECT_TRUE(ss_residue.CA());
    EXPECT_TRUE(ss_residue.C());
    EXPECT_TRUE(ss_residue.O());
}

TEST_F(SSResidueTest, FromNonProteinResidue)
{
    SSResidue ss_residue(non_protein_residue);

    EXPECT_FALSE(ss_residue.is_amino_acid());
    EXPECT_EQ(ss_residue.secondary_structure(), mol::Unknown);
    EXPECT_FALSE(ss_residue.is_proline());
    EXPECT_FALSE(ss_residue.is_chain_break());
    EXPECT_EQ(ss_residue.chain(), "C");
    EXPECT_FALSE(ss_residue.N());
    EXPECT_FALSE(ss_residue.CA());
    EXPECT_FALSE(ss_residue.C());
    EXPECT_FALSE(ss_residue.O());
}

TEST_F(SSResidueTest, DefaultFrame)
{
    SSResidue ss_protein_residue(protein_residue);
    SSResidue ss_proline_residue(proline_residue);
    SSResidue ss_non_protein_residue(non_protein_residue);

    EXPECT_EQ(ss_protein_residue.N().frame(), 0);
    EXPECT_EQ(ss_protein_residue.CA().frame(), 0);
    EXPECT_EQ(ss_protein_residue.C().frame(), 0);
    EXPECT_EQ(ss_protein_residue.O().frame(), 0);

    EXPECT_EQ(ss_proline_residue.N().frame(), 0);
    EXPECT_EQ(ss_proline_residue.CA().frame(), 0);
    EXPECT_EQ(ss_proline_residue.C().frame(), 0);
    EXPECT_EQ(ss_proline_residue.O().frame(), 0);

    EXPECT_FALSE(ss_non_protein_residue.N().frame());
    EXPECT_FALSE(ss_non_protein_residue.CA().frame());
    EXPECT_FALSE(ss_non_protein_residue.C().frame());
    EXPECT_FALSE(ss_non_protein_residue.O().frame());
}

TEST_F(SSResidueTest, SetValidFrame)
{
    SSResidue ss_protein_residue(protein_residue);
    SSResidue ss_proline_residue(proline_residue);
    SSResidue ss_non_protein_residue(non_protein_residue);

    ss_protein_residue.set_frame(1);
    ss_proline_residue.set_frame(1);
    ss_non_protein_residue.set_frame(1);

    EXPECT_EQ(ss_protein_residue.N().frame(), 1);
    EXPECT_EQ(ss_protein_residue.CA().frame(), 1);
    EXPECT_EQ(ss_protein_residue.C().frame(), 1);
    EXPECT_EQ(ss_protein_residue.O().frame(), 1);

    EXPECT_EQ(ss_protein_residue.N().frame(), 1);
    EXPECT_EQ(ss_protein_residue.CA().frame(), 1);
    EXPECT_EQ(ss_protein_residue.C().frame(), 1);
    EXPECT_EQ(ss_protein_residue.O().frame(), 1);

    EXPECT_FALSE(ss_non_protein_residue.N().frame());
    EXPECT_FALSE(ss_non_protein_residue.CA().frame());
    EXPECT_FALSE(ss_non_protein_residue.C().frame());
    EXPECT_FALSE(ss_non_protein_residue.O().frame());
}

TEST_F(SSResidueTest, SetInvalidFrame)
{
    SSResidue ss_protein_residue(protein_residue);
    SSResidue ss_proline_residue(proline_residue);
    SSResidue ss_non_protein_residue(non_protein_residue);

    EXPECT_THROW(ss_protein_residue.set_frame(2), MolError);
    EXPECT_THROW(ss_proline_residue.set_frame(2), MolError);
    EXPECT_NO_THROW(ss_non_protein_residue.set_frame(2));

    EXPECT_FALSE(ss_non_protein_residue.N().frame());
    EXPECT_FALSE(ss_non_protein_residue.CA().frame());
    EXPECT_FALSE(ss_non_protein_residue.C().frame());
    EXPECT_FALSE(ss_non_protein_residue.O().frame());
}

TEST_F(SSResidueTest, SetSecondaryStructure)
{
    SSResidue ss_protein_residue(protein_residue);
    SSResidue ss_proline_residue(proline_residue);
    SSResidue ss_non_protein_residue(non_protein_residue);

    EXPECT_EQ(ss_protein_residue.secondary_structure(), mol::Unknown);
    EXPECT_EQ(ss_proline_residue.secondary_structure(), mol::Unknown);
    EXPECT_EQ(ss_non_protein_residue.secondary_structure(), mol::Unknown);

    ss_protein_residue.set_secondary_structure(mol::Helix3);
    ss_proline_residue.set_secondary_structure(mol::Helix3);
    ss_non_protein_residue.set_secondary_structure(mol::Helix3);

    EXPECT_EQ(ss_protein_residue.secondary_structure(), mol::Helix3);
    EXPECT_EQ(ss_proline_residue.secondary_structure(), mol::Helix3);
    EXPECT_EQ(ss_non_protein_residue.secondary_structure(), mol::Helix3);
}

TEST_F(SSResidueTest, SetChainBreak)
{
    SSResidue ss_protein_residue(protein_residue);
    SSResidue ss_proline_residue(proline_residue);
    SSResidue ss_non_protein_residue(non_protein_residue);

    EXPECT_FALSE(ss_protein_residue.is_chain_break());
    EXPECT_FALSE(ss_proline_residue.is_chain_break());
    EXPECT_FALSE(ss_non_protein_residue.is_chain_break());

    ss_protein_residue.set_chain_break(true);
    ss_proline_residue.set_chain_break(true);
    ss_non_protein_residue.set_chain_break(true);

    EXPECT_TRUE(ss_protein_residue.is_chain_break());
    EXPECT_TRUE(ss_proline_residue.is_chain_break());
    EXPECT_TRUE(ss_non_protein_residue.is_chain_break());
}

class DSSPTest : public ::testing::Test
{
public:
    static MolSystem const& all_structures()
    {
        static MolSystem mol = load_mol("s_structure.pdb");
        return mol;
    }

    static MolSystem const& multiple_molecules()
    {
        static MolSystem mol = load_mol("4lad.pdb");
        return mol;
    }

    constexpr std::array<mol::SecondaryStructure, 585> ss_all_structures();
    constexpr std::array<mol::SecondaryStructure, 266> ss_multiple_molecules();

private:
    static MolSystem load_mol(std::string const& file_name)
    {
        MolSystem mol(file_name);
        mol.add_trajectory(file_name);
        return mol;
    }
};

TEST_F(DSSPTest, ProteinAndNonProtein)
{
    DSSP dssp(multiple_molecules());
    dssp.run(0);
    EXPECT_THAT(dssp.secondary_structures(), ElementsAreArray(ss_multiple_molecules()));
    EXPECT_THAT(dssp.residues(), Pointwise(Prop(&SSResidue::secondary_structure), ss_multiple_molecules()));
}

TEST_F(DSSPTest, MultipleRuns)
{
    DSSP dssp(all_structures());
    dssp.run(0);
    dssp.run(0);
    EXPECT_THAT(dssp.secondary_structures(), ElementsAreArray(ss_all_structures()));
    EXPECT_THAT(dssp.residues(), Pointwise(Prop(&SSResidue::secondary_structure), ss_all_structures()));
}

TEST_F(DSSPTest, AllStructures)
{
    DSSP dssp(all_structures());
    dssp.run(0);
    EXPECT_THAT(dssp.secondary_structures(), ElementsAreArray(ss_all_structures()));
    EXPECT_THAT(dssp.residues(), Pointwise(Prop(&SSResidue::secondary_structure), ss_all_structures()));
}

TEST_F(DSSPTest, InvalidFrame)
{
    DSSP dssp(all_structures());
    EXPECT_THROW(dssp.run(1), MolError);
}

TEST_F(DSSPTest, ProlineDetection)
{
    DSSP dssp(all_structures());
    dssp.run(0);
    for (size_t const index : std::vector<size_t>{44, 56, 135, 210, 215, 234, 254, 278, 295, 298, 335, 339, 357, 372, 374, 440, 460, 469, 508, 509, 516, 532, 543, 549, 553, 556, 557, 568, 572, 584})
    {
        EXPECT_TRUE(dssp.residues()[index].is_proline()) << index;
    }
}

TEST_F(DSSPTest, ChainBreakDetection)
{
    DSSP dssp(all_structures());
    dssp.run(0);
    auto res = dssp.residues();
    for (size_t const index : std::vector<size_t>{179, 491})
    {
        EXPECT_TRUE(dssp.residues()[index].is_chain_break()) << index;
    }
}

constexpr std::array<mol::SecondaryStructure, 585> DSSPTest::ss_all_structures()
{
    std::array<mol::SecondaryStructure, 585> structure{
        mol::Loop, mol::Loop, mol::Bend, mol::Bend, mol::Loop, mol::Loop, mol::Loop, mol::Loop, mol::Helix3, mol::Helix3, mol::Helix3, mol::Helix3, mol::Helix3, mol::Helix3, mol::Helix3, mol::Loop, mol::Loop, mol::Loop, mol::Loop, mol::Loop, mol::Loop, mol::Loop, mol::Bend, mol::Bend, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix5, mol::Helix5, mol::Helix5, mol::Helix5, mol::Helix5, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Turn, mol::Turn, mol::Loop, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Bend, mol::Turn, mol::Turn, mol::Loop, mol::Loop, mol::Loop, mol::Loop, mol::Loop, mol::Loop, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Turn, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Bend, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix5, mol::Helix5, mol::Helix5, mol::Helix5, mol::Helix5, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Turn, mol::Turn, mol::Loop, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Loop, mol::Loop, mol::Loop, mol::Loop, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Bend, mol::Loop, mol::Turn, mol::Turn, mol::Loop, mol::Helix3, mol::Helix3, mol::Helix3, mol::Helix3, mol::Bend, mol::Loop, mol::Bend, mol::Bend, mol::Bend, mol::Loop, mol::Loop, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Bend, mol::Loop, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Turn, mol::Loop, mol::Bridge, mol::Turn, mol::Turn, mol::Turn, mol::Turn, mol::Bridge, mol::Loop, mol::Loop, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Turn, mol::Bend, mol::Loop, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Bend, mol::Loop, mol::Turn, mol::Turn, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Turn, mol::Bend, mol::Loop, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Loop, mol::Loop, mol::Bend, mol::Bend, mol::Loop, mol::Loop, mol::Loop, mol::Loop, mol::Loop, mol::Loop, mol::Loop, mol::Loop, mol::Loop, mol::Bend, mol::Loop, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Bend, mol::Loop, mol::Turn, mol::Turn, mol::Loop, mol::Bend, mol::Helix3, mol::Helix3, mol::Helix3, mol::Bend, mol::Loop, mol::Loop, mol::Turn, mol::Turn, mol::Loop, mol::Loop, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Bend, mol::Loop, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Turn, mol::Loop, mol::Loop, mol::Turn, mol::Turn, mol::Bend, mol::Loop, mol::Loop, mol::Loop, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Turn, mol::Turn, mol::Bend, mol::Loop, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Loop, mol::Loop, mol::Turn, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Turn, mol::Turn, mol::Bend, mol::Loop, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Loop, mol::Loop, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Bend, mol::Loop, mol::Loop, mol::Turn, mol::Turn, mol::Strand, mol::Strand, mol::Strand, mol::Strand, mol::Strand, mol::Bend, mol::Bend, mol::Turn, mol::Turn, mol::Loop, mol::Turn, mol::Turn, mol::Strand, mol::Strand, mol::Strand, mol::Strand, mol::Strand, mol::Strand, mol::Strand, mol::Loop, mol::Loop, mol::Turn, mol::Turn, mol::Bend, mol::Turn, mol::Turn, mol::Turn, mol::Turn, mol::Loop, mol::Strand, mol::Strand, mol::Strand, mol::Strand, mol::Strand, mol::Strand, mol::Strand, mol::Loop, mol::Loop, mol::Turn, mol::Turn, mol::Turn, mol::Turn, mol::Bend, mol::Bend, mol::Loop, mol::Loop, mol::Strand, mol::Strand, mol::Strand, mol::Strand, mol::Loop, mol::Bend, mol::Loop, mol::Loop, mol::Loop, mol::Loop, mol::Turn, mol::Turn, mol::Bridge, mol::Loop, mol::Turn, mol::Turn, mol::Loop, mol::Bridge, mol::Bridge, mol::Loop, mol::Loop, mol::Helix3, mol::Helix3, mol::Helix3, mol::Bend, mol::Loop, mol::Loop,
    };

    return structure;
}

constexpr std::array<mol::SecondaryStructure, 266> DSSPTest::ss_multiple_molecules()
{
    std::array<mol::SecondaryStructure, 266> structure{
        mol::Loop, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Bend, mol::Loop, mol::Loop, mol::Turn, mol::Turn, mol::Strand, mol::Strand, mol::Strand, mol::Strand, mol::Strand, mol::Bend, mol::Bend, mol::Turn, mol::Turn, mol::Loop, mol::Turn, mol::Turn, mol::Strand, mol::Strand, mol::Strand, mol::Strand, mol::Strand, mol::Strand, mol::Strand, mol::Loop, mol::Loop, mol::Turn, mol::Turn, mol::Bend, mol::Turn, mol::Turn, mol::Turn, mol::Turn, mol::Loop, mol::Strand, mol::Strand, mol::Strand, mol::Strand, mol::Strand, mol::Strand, mol::Strand, mol::Loop, mol::Loop, mol::Turn, mol::Turn, mol::Turn, mol::Turn, mol::Bend, mol::Bend, mol::Loop, mol::Loop, mol::Strand, mol::Strand, mol::Strand, mol::Strand, mol::Loop, mol::Bend, mol::Loop, mol::Loop, mol::Loop, mol::Loop, mol::Turn, mol::Turn, mol::Bridge, mol::Loop, mol::Turn, mol::Turn, mol::Loop, mol::Bridge, mol::Bridge, mol::Loop, mol::Loop, mol::Helix3, mol::Helix3, mol::Helix3, mol::Bend, mol::Loop, mol::Loop, mol::Loop, mol::Loop, mol::Loop, mol::Loop, mol::Loop, mol::Loop, mol::Turn, mol::Turn, mol::Loop, mol::Loop, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Bend, mol::Loop, mol::Loop, mol::Bend, mol::Bend, mol::Loop, mol::Loop, mol::Bend, mol::Loop, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Turn, mol::Turn, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Turn, mol::Turn, mol::Loop, mol::Loop, mol::Loop, mol::Loop, mol::Bend, mol::Bend, mol::Bend, mol::Loop, mol::Bend, mol::Loop, mol::Loop, mol::Loop, mol::Loop, mol::Strand, mol::Strand, mol::Strand, mol::Loop, mol::Turn, mol::Turn, mol::Bend, mol::Loop, mol::Strand, mol::Strand, mol::Strand, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Turn, mol::Turn, mol::Loop, mol::Bend, mol::Bend, mol::Loop, mol::Bend, mol::Bend, mol::Bend, mol::Loop, mol::Loop, mol::Loop, mol::Loop, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Turn, mol::Loop, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown,
    };

    return structure;
}

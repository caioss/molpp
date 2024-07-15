#include "guessers/ResidueBondGuesser.hpp"
#include <molpp/internal/MolData.hpp>
#include <molpp/Atom.hpp>
#include <molpp/Residue.hpp>
#include <molpp/AtomSel.hpp>
#include <molpp/ResidueSel.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <array>

using namespace mol;
using namespace mol::internal;
using namespace testing;

//! Test fixture for ResidueBondGuesser.
class ResidueBondGuesserTest : public testing::Test
{
public:
    ResidueBondGuesserTest()
    : data{create_moldata()}
    , atoms{std::to_array<mol::index_t>({0, 1, 2, 3, 4, 5, 6, 7, 8}), data}
    , residues{std::to_array<mol::index_t>({0, 1, 2}), data}
    {
        guesser.apply(residues);
    }

    MolData data;
    AtomSel atoms;
    ResidueSel residues;
    ResidueBondGuesser guesser;

private:
    MolData create_moldata()
    {
        MolData data(9);
        ResidueData& res_data = data.residues();
        res_data.resize(3);
        data.topology().link_categories(EntityCategory::Residue, EntityCategory::Atom);

        // GLY
        Atom glycine_N(0, std::nullopt, data);
        glycine_N.set_name("N");
        Atom(1, std::nullopt, data).set_name("CA");
        Atom(2, std::nullopt, data).set_name("C");
        Residue glycine(0, std::nullopt, data);
        glycine.set_name("GLY");
        glycine.add_atom(0);
        glycine.add_atom(1);
        glycine.add_atom(2);

        // Add existing bond
        auto bond = glycine_N.add_bond(1);
        bond->set_aromatic(true);
        bond->set_guessed(false);
        bond->set_order(2);
        bond->set_guessed_order(false);

        // MET
        Atom(3, std::nullopt, data).set_name("N");
        Atom(4, std::nullopt, data).set_name("CG");
        Atom(5, std::nullopt, data).set_name("SD");
        Atom(6, std::nullopt, data).set_name("UNK");
        Residue metionine(1, std::nullopt, data);
        metionine.set_name("MET");
        metionine.add_atom(3);
        metionine.add_atom(4);
        metionine.add_atom(5);
        metionine.add_atom(6);

        // Modified metionine
        Atom(7, std::nullopt, data).set_name("CG");
        Atom(8, std::nullopt, data).set_name("SD");
        Residue modified(2, std::nullopt, data);
        modified.set_name("MMT");
        modified.add_atom(7);
        modified.add_atom(8);

        return data;
    }
};

// Guessed bonds should be added.
TEST_F(ResidueBondGuesserTest, guess_bonds)
{
    EXPECT_THAT(atoms[1].bond(2), NotNull());
    EXPECT_THAT(atoms[4].bond(5), NotNull());
}

// Existing bonds should not be modified.
TEST_F(ResidueBondGuesserTest, existing_bonds)
{
    auto bond = atoms[0].bond(1);

    ASSERT_THAT(bond, NotNull());
    EXPECT_EQ(bond->order(), 2);
    EXPECT_TRUE(bond->aromatic());
    EXPECT_FALSE(bond->guessed());
    EXPECT_FALSE(bond->guessed_order());
}

// Bonds should not be added to peptide bonds.
TEST_F(ResidueBondGuesserTest, peptide_bonds)
{
    EXPECT_THAT(atoms[2].bond(3), IsNull());
}

// Unknown atoms of known residues should not be bonded.
TEST_F(ResidueBondGuesserTest, unknown_atoms)
{
    EXPECT_THAT(atoms[6].bonds(), IsEmpty());
}

// Known atoms of unknown residues should not be bonded.
TEST_F(ResidueBondGuesserTest, unknown_residues)
{
    EXPECT_THAT(atoms[7].bonds(), IsEmpty());
    EXPECT_THAT(atoms[8].bonds(), IsEmpty());
}

// Atoms with missing bond partners should not be bonded.
TEST_F(ResidueBondGuesserTest, missing_bond_partners)
{
    EXPECT_THAT(atoms[3].bonds(), IsEmpty());
}

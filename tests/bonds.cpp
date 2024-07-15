#include <molpp/AtomSel.hpp>
#include <molpp/ResidueSel.hpp>
#include <molpp/internal/MolData.hpp>
#include "guessers/ResidueBondGuesser.hpp"
#include "guessers/AtomBondGuesser.hpp"
#include "files.hpp"
#include "matchers.hpp"
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <vector>

using namespace mol;
using namespace mol::internal;
using namespace testing;

TEST(Bonds, BondData) {
    BondData bond_data(5);

    EXPECT_TRUE(bond_data.incomplete());
    bond_data.set_incomplete(false);
    EXPECT_FALSE(bond_data.incomplete());

    EXPECT_THAT(bond_data.add_bond(1, 2), NotNull());
    EXPECT_EQ(bond_data.add_bond(2, 1), bond_data.bond(1, 2));
    EXPECT_THAT(bond_data.add_bond(3, 1), NotNull());
    EXPECT_THAT(bond_data.add_bond(4, 3), NotNull());
    EXPECT_THAT(bond_data.add_bond(1, 1), IsNull());

    EXPECT_THAT(bond_data.bonded(0), UnorderedElementsAre());
    EXPECT_THAT(bond_data.bonded(1), UnorderedElementsAre(1, 2, 3));
    EXPECT_THAT(bond_data.bonded(2), UnorderedElementsAre(1, 2));
    EXPECT_THAT(bond_data.bonded(3), UnorderedElementsAre(1, 3, 4));
    EXPECT_THAT(bond_data.bonded(4), ElementsAre(3, 4));

    index_t indices[3] = {1, 0, 4};
    auto bonded = bond_data.bonds(indices, indices + 3);
    EXPECT_EQ(bonded.size(), 3);
    std::vector<index_t> bond_indices;
    bond_indices.reserve(6);
    for (auto b : bonded)
    {
        bond_indices.push_back(b->atom1());
        bond_indices.push_back(b->atom2());
    }
    EXPECT_THAT(bond_indices, UnorderedElementsAre(1, 2, 3, 1, 3, 4));
    EXPECT_THAT(bond_data.bonded(indices, indices + 3), UnorderedElementsAre(1, 2, 3, 4));

    bonded = bond_data.bonds(1);
    EXPECT_EQ(bonded.size(), 2);
    bond_indices.clear();
    bond_indices.reserve(4);
    for (auto b : bonded)
    {
        bond_indices.push_back(b->atom1());
        bond_indices.push_back(b->atom2());
    }
    EXPECT_THAT(bond_indices, UnorderedElementsAre(1, 1, 2, 3));

    EXPECT_THAT(bond_data.bond(1, 5), IsNull());
    auto bond = bond_data.bond(2, 1);
    EXPECT_THAT(bond, NotNull());
    EXPECT_EQ(bond, bond_data.bond(1, 2));
    EXPECT_EQ(bond->atom1(), 1);
    EXPECT_EQ(bond->atom2(), 2);
    EXPECT_EQ(bond->order(), 0);
    EXPECT_FALSE(bond->aromatic());
    EXPECT_TRUE(bond->guessed());
    EXPECT_TRUE(bond->guessed_order());

    EXPECT_EQ(bond_data.size(), 3);
    bond_data.clear();
    EXPECT_EQ(bond_data.size(), 0);
}

TEST(Bonds, Bond) {
    Bond bond(2, 1);

    EXPECT_EQ(bond.order(), 0);
    bond.set_order(1);
    EXPECT_EQ(bond.order(), 1);

    EXPECT_TRUE(bond.guessed());
    bond.set_guessed(false);
    EXPECT_FALSE(bond.guessed());

    EXPECT_TRUE(bond.guessed_order());
    bond.set_guessed_order(false);
    EXPECT_FALSE(bond.guessed_order());

    EXPECT_FALSE(bond.aromatic());
    bond.set_aromatic(true);
    EXPECT_TRUE(bond.aromatic());

    EXPECT_EQ(bond.atom1(), 1);
    EXPECT_EQ(bond.atom2(), 2);
    bond = Bond(1, 2);
    EXPECT_EQ(bond.atom1(), 1);
    EXPECT_EQ(bond.atom2(), 2);
}

TEST(Bonds, Guessers) {
    std::shared_ptr<MolReader> reader = MolReader::from_file_ext(".pdb");
    ASSERT_THAT(reader, NotNull());

    /*
     * Atom-based guesser
     */
    auto atom_data = reader->read_topology("4lad.pdb");
    ASSERT_THAT(atom_data, NotNull());
    atom_data->bonds().clear();
    reader->read_trajectory("4lad.pdb", *atom_data);

    AtomSel atoms_sel(*atom_data);
    AtomBondGuesser atom_guesser;
    atom_guesser.apply(atoms_sel);

    // Peptide bond
    auto atom_bond = atoms_sel[705].bond(711); // ILE90A-C-SER91A-N
    ASSERT_THAT(atom_bond, NotNull());
    EXPECT_EQ(atom_bond->order(), 1);
    EXPECT_FALSE(atom_bond->aromatic());
    EXPECT_TRUE(atom_bond->guessed());
    EXPECT_TRUE(atom_bond->guessed_order());

    // Metals
    atom_bond = atoms_sel[1791].bond(1407); // HIS361B-ND1-ZN701B
    ASSERT_THAT(atom_bond, NotNull());
    EXPECT_EQ(atom_bond->order(), 1);
    EXPECT_FALSE(atom_bond->aromatic());
    EXPECT_TRUE(atom_bond->guessed());
    EXPECT_TRUE(atom_bond->guessed_order());

    atom_bond = atoms_sel[1791].bond(1430); // CYS364B-SG-ZN701B
    ASSERT_THAT(atom_bond, NotNull());
    EXPECT_EQ(atom_bond->order(), 1);
    EXPECT_FALSE(atom_bond->aromatic());
    EXPECT_TRUE(atom_bond->guessed());
    EXPECT_TRUE(atom_bond->guessed_order());

    // Extra molecules
    atom_bond = atoms_sel[1798].bond(1794); // OXL703B
    ASSERT_THAT(atom_bond, NotNull());
    EXPECT_EQ(atom_bond->order(), 1);
    EXPECT_FALSE(atom_bond->aromatic());
    EXPECT_TRUE(atom_bond->guessed());
    EXPECT_TRUE(atom_bond->guessed_order());
    atom_bond = atoms_sel[1796].bond(1794); // OXL703B
    ASSERT_THAT(atom_bond, NotNull());
    EXPECT_EQ(atom_bond->order(), 1);
    EXPECT_FALSE(atom_bond->aromatic());
    EXPECT_TRUE(atom_bond->guessed());
    EXPECT_TRUE(atom_bond->guessed_order());
    atom_bond = atoms_sel[1793].bond(1794); // OXL703B
    ASSERT_THAT(atom_bond, NotNull());
    EXPECT_EQ(atom_bond->order(), 1);
    EXPECT_FALSE(atom_bond->aromatic());
    EXPECT_TRUE(atom_bond->guessed());
    EXPECT_TRUE(atom_bond->guessed_order());
    atom_bond = atoms_sel[1793].bond(1797); // OXL703B
    ASSERT_THAT(atom_bond, NotNull());
    EXPECT_EQ(atom_bond->order(), 1);
    EXPECT_FALSE(atom_bond->aromatic());
    EXPECT_TRUE(atom_bond->guessed());
    EXPECT_TRUE(atom_bond->guessed_order());
    atom_bond = atoms_sel[1793].bond(1795); // OXL703B
    ASSERT_THAT(atom_bond, NotNull());
    EXPECT_EQ(atom_bond->order(), 1);
    EXPECT_FALSE(atom_bond->aromatic());
    EXPECT_TRUE(atom_bond->guessed());
    EXPECT_TRUE(atom_bond->guessed_order());

    // Water
    atom_bond = atoms_sel[978].bond(1830); // GLY135A-N-HOH226
    EXPECT_THAT(atom_bond, IsNull());
}

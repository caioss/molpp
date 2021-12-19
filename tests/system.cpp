#include <molpp/Atom.hpp>
#include <molpp/MolSystem.hpp>
#include <molpp/MolError.hpp>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace testing;

TEST(System, MolSystem) {
    EXPECT_THROW(MolSystem(""), MolError);
    EXPECT_THROW(MolSystem("file.unk"), MolError);
    EXPECT_THROW(MolSystem("no_file.pdb"), MolError);
    MolSystem mol("traj.pdb");

    EXPECT_THROW(mol.add_trajectory("traj.unk"), MolError);
    EXPECT_THROW(mol.add_trajectory("tiny.pdb"), MolError);
    EXPECT_THROW(mol.add_trajectory("no_file.pdb"), MolError);
    EXPECT_NO_THROW(mol.add_trajectory("traj.pdb"));

    AtomSel all_sel{mol.atoms()};
    EXPECT_EQ(all_sel.size(), 2);

    AtomSel index_sel{mol.select({0})};
    EXPECT_EQ(index_sel.size(), 1);
}

TEST(System, Bonds) {
    MolSystem mol("4lad.pdb");
    mol.add_trajectory("4lad.pdb");
    AtomSel atoms{mol.atoms()};

    ASSERT_EQ(atoms.size(), 1844);
    ASSERT_TRUE(atoms.frame());

    // Pre-condition
    auto bond = atoms[650].bond(648); // TYR80A-CZ-CE1
    EXPECT_THAT(bond, IsNull());
    bond = atoms[1108].bond(1101); // PHE151A-N-GLN150A-C
    EXPECT_THAT(bond, IsNull());
    bond = atoms[1798].bond(1794); // OXL703B
    EXPECT_THAT(bond, NotNull());
    bond = atoms[1791].bond(1407); // HIS361B-ND1-ZN701B
    EXPECT_THAT(bond, NotNull());

    // Reset bonds
    mol.reset_bonds();
    bond = atoms[650].bond(648); // TYR80A-CZ-CE1
    EXPECT_THAT(bond, IsNull());
    bond = atoms[1108].bond(1101); // PHE151A-N-GLN150A-C
    EXPECT_THAT(bond, IsNull());
    bond = atoms[1798].bond(1794); // OXL703B
    EXPECT_THAT(bond, IsNull());
    bond = atoms[1791].bond(1407); // HIS361B-ND1-ZN701B
    EXPECT_THAT(bond, IsNull());

    // Guess bonds
    mol.guess_bonds();
    bond = atoms[650].bond(648); // TYR80A-CZ-CE1
    EXPECT_THAT(bond, NotNull());
    bond = atoms[1108].bond(1101); // PHE151A-N-GLN150A-C
    EXPECT_THAT(bond, NotNull());
    bond = atoms[1798].bond(1794); // OXL703B
    EXPECT_THAT(bond, NotNull());
    bond = atoms[1791].bond(1407); // HIS361B-ND1-ZN701B
    EXPECT_THAT(bond, NotNull());
}

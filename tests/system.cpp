#include <molpp/Atom.hpp>
#include <molpp/MolError.hpp>
#include <molpp/MolSystem.hpp>
#include <molpp/AtomSelector.hpp>
#include <molpp/internal/MolData.hpp>
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

}

TEST(System, Selection) {
    MolSystem mol("4lad.pdb");
    mol.add_trajectory("4lad.pdb");
    mol.add_trajectory("4lad.pdb");

    // All atoms selections
    AtomSel all_sel{mol.atoms()};
    EXPECT_EQ(all_sel.size(), 1844);
    EXPECT_FALSE(all_sel.frame());
    all_sel = mol.select("all");
    EXPECT_EQ(all_sel.size(), 1844);
    EXPECT_FALSE(all_sel.frame());

    // All atoms with frame
    EXPECT_EQ(mol.atoms(1).frame(), 1);
    EXPECT_EQ(mol.select("all", 1).frame(), 1);

    // From indices with frame
    AtomSel index_sel{mol.select(std::vector<index_t>{0, 2}, 1)};
    EXPECT_THAT(index_sel.indices(), ElementsAre(0, 2));
    EXPECT_EQ(index_sel.frame(), 1);
    EXPECT_EQ(mol.select(std::vector<index_t>{0, 2}, 1).frame(), 1);

    // From indices without frame
    index_sel = mol.select(std::vector<index_t>{0, 2});
    EXPECT_THAT(index_sel.indices(), ElementsAre(0, 2));
    EXPECT_FALSE(index_sel.frame());
    EXPECT_EQ(mol.select(std::vector<index_t>{0, 2}, 1).frame(), 1);

    // From strings with frame
    AtomSel water = mol.select("resid 203:205 or resid 900:910", 1);
    EXPECT_THAT(water.indices(), ElementsAre(1807, 1808, 1809));
    EXPECT_EQ(water.frame(), 1);

    // From strings without frame
    water = mol.select("resid 203:205 or resid 900:910");
    EXPECT_THAT(water.indices(), ElementsAre(1807, 1808, 1809));
    EXPECT_FALSE(water.frame());

    // Selector
    AtomSelector selector = mol.selector("resid 203:205 or resid 900:910");
    water = selector.apply(1);
    EXPECT_THAT(water.indices(), ElementsAre(1807, 1808, 1809));
    EXPECT_EQ(water.frame(), 1);
}

TEST(System, Bonds) {
    MolSystem mol("4lad.pdb");
    mol.add_trajectory("4lad.pdb");
    AtomSel atoms{mol.atoms(0)};

    ASSERT_EQ(atoms.size(), 1844);
    ASSERT_EQ(atoms.frame(), 0);

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
    mol.guess_bonds(0);
    bond = atoms[650].bond(648); // TYR80A-CZ-CE1
    EXPECT_THAT(bond, NotNull());
    bond = atoms[1108].bond(1101); // PHE151A-N-GLN150A-C
    EXPECT_THAT(bond, NotNull());
    bond = atoms[1798].bond(1794); // OXL703B
    EXPECT_THAT(bond, NotNull());
    bond = atoms[1791].bond(1407); // HIS361B-ND1-ZN701B
    EXPECT_THAT(bond, NotNull());
}

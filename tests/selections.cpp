#include "files.hpp"
#include "matchers.hpp"
#include <molpp/Atom.hpp>
#include <molpp/AtomSel.hpp>
#include <molpp/MolError.hpp>
#include <molpp/internal/BaseSel.hpp>
#include <molpp/internal/Sel.hpp>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace mol::internal;
using namespace testing;

TEST(Selections, BaseSel) {
    PDBFiles pdb;
    pdb.check();

    /*
     * Construction
     */
    BaseSel all_sel(pdb.tiny->size(), pdb.tiny);
    EXPECT_EQ(all_sel.size(), pdb.tiny->size());
    EXPECT_FALSE(all_sel.frame());
    EXPECT_THAT(all_sel.indices(), ElementsAre(0, 1, 2, 3, 4, 5));
    EXPECT_THAT(all_sel.selected(), ElementsAre(true, true, true, true, true, true));
    for (size_t i = 0; i < 6; ++i)
    {
        EXPECT_TRUE(all_sel.contains(i)) << "index " << i;
    }
    EXPECT_FALSE(all_sel.contains(6));

    // Constructors accepting indexes
    std::vector<size_t> indices{4, 1, 3};
    BaseSel some_sel(pdb.tiny->size(), indices, pdb.tiny);
    BaseSel rvalue_sel(pdb.tiny->size(), {4, 1, 3}, pdb.tiny);
    EXPECT_EQ(some_sel.size(), indices.size());
    EXPECT_EQ(rvalue_sel.size(), indices.size());
    EXPECT_FALSE(some_sel.frame());
    EXPECT_FALSE(rvalue_sel.frame());
    EXPECT_THAT(some_sel.indices(), ElementsAre(1, 3, 4));
    EXPECT_THAT(rvalue_sel.indices(), ElementsAre(1, 3, 4));
    EXPECT_THAT(some_sel.selected(), ElementsAre(false, true, false, true, true, false));
    EXPECT_THAT(rvalue_sel.selected(), ElementsAre(false, true, false, true, true, false));
    for (size_t i : {1, 3, 4})
    {
        EXPECT_TRUE(some_sel.contains(i)) << "Index " << i;
        EXPECT_TRUE(rvalue_sel.contains(i)) << "Index " << i;
    }
    for (size_t i : {0, 2, 5})
    {
        EXPECT_FALSE(some_sel.contains(i)) << "Index " << i;
        EXPECT_FALSE(rvalue_sel.contains(i)) << "Index " << i;
    }

    /*
     * Trajectory
     */
    BaseSel traj_sel(pdb.traj->size(), pdb.traj);
    ASSERT_TRUE(traj_sel.frame());
    EXPECT_EQ(traj_sel.frame(), 0);
    traj_sel.set_frame(3);
    ASSERT_TRUE(traj_sel.frame());
    EXPECT_EQ(traj_sel.frame(), 3);
    EXPECT_THROW(traj_sel.set_frame(4), MolError);
    ASSERT_TRUE(traj_sel.frame());
    EXPECT_EQ(traj_sel.frame(), 3);
}

TEST(Selections, Sel) {
    PDBFiles pdb;
    pdb.check();

    /*
     * Construction
     */
    Sel<Atom, AtomSel> all_sel(pdb.tiny);
    EXPECT_FALSE(all_sel.frame());
    EXPECT_EQ(all_sel.size(), pdb.tiny->size());
    EXPECT_THAT(all_sel.indices(), ElementsAre(0, 1, 2, 3, 4, 5));
    EXPECT_THAT(all_sel.selected(), ElementsAre(true, true, true, true, true, true));
    for (size_t i = 0; i < 6; ++i)
    {
        EXPECT_TRUE(all_sel.contains(i)) << "index " << i;
    }
    EXPECT_FALSE(all_sel.contains(6));

    // Constructors accepting indexes
    std::vector<size_t> indices{4, 1, 3};
    Sel<Atom, AtomSel> some_sel(indices, pdb.tiny);
    Sel<Atom, AtomSel> rvalue_sel({4, 1, 3}, pdb.tiny);
    EXPECT_EQ(some_sel.size(), indices.size());
    EXPECT_EQ(rvalue_sel.size(), indices.size());
    EXPECT_FALSE(some_sel.frame());
    EXPECT_FALSE(rvalue_sel.frame());
    EXPECT_THAT(some_sel.indices(), ElementsAre(1, 3, 4));
    EXPECT_THAT(rvalue_sel.indices(), ElementsAre(1, 3, 4));
    EXPECT_THAT(some_sel.selected(), ElementsAre(false, true, false, true, true, false));
    EXPECT_THAT(rvalue_sel.selected(), ElementsAre(false, true, false, true, true, false));
    for (size_t i : {1, 3, 4})
    {
        EXPECT_TRUE(some_sel.contains(i)) << "Index " << i;
        EXPECT_TRUE(rvalue_sel.contains(i)) << "Index " << i;
    }
    for (size_t i : {0, 2, 5})
    {
        EXPECT_FALSE(some_sel.contains(i)) << "Index " << i;
        EXPECT_FALSE(rvalue_sel.contains(i)) << "Index " << i;
    }

    /*
     * Atoms getters
     */
    EXPECT_EQ(some_sel[0].index(), 1);
    EXPECT_FALSE(some_sel[0].frame());
    EXPECT_EQ(some_sel.at(0).index(), 1);
    EXPECT_FALSE(some_sel.at(0).frame());
    EXPECT_EQ(rvalue_sel[0].index(), 1);
    EXPECT_EQ(some_sel[1].index(), 3);
    EXPECT_FALSE(some_sel[1].frame());
    EXPECT_EQ(some_sel.at(1).index(), 3);
    EXPECT_FALSE(some_sel.at(1).frame());
    EXPECT_EQ(rvalue_sel[1].index(), 3);
    EXPECT_EQ(some_sel[2].index(), 4);
    EXPECT_FALSE(some_sel[2].frame());
    EXPECT_EQ(some_sel.at(2).index(), 4);
    EXPECT_FALSE(some_sel.at(2).frame());
    EXPECT_EQ(rvalue_sel[2].index(), 4);
    EXPECT_THROW(some_sel.at(3).index(), MolError);
    for (size_t i : {1, 3, 4})
    {
        EXPECT_TRUE(some_sel.contains(i)) << "Index " << i;
        EXPECT_TRUE(rvalue_sel.contains(i)) << "Index " << i;
    }
    for (size_t i : {0, 2, 5})
    {
        EXPECT_FALSE(some_sel.contains(i)) << "Index " << i;
        EXPECT_FALSE(rvalue_sel.contains(i)) << "Index " << i;
    }

    /*
     * Trajectory
     */
    EXPECT_THROW(all_sel.coords(), MolError);
    Sel<Atom, AtomSel> traj_sel(pdb.traj);
    ASSERT_TRUE(traj_sel.frame());
    EXPECT_EQ(traj_sel.frame(), 0);
    EXPECT_THAT(traj_sel.coords().reshaped(), ElementsAre(1, -1, 0, 2, -2, 1));
    traj_sel.set_frame(3);
    ASSERT_TRUE(traj_sel.frame());
    EXPECT_EQ(traj_sel.frame(), 3);
    EXPECT_THROW(traj_sel.set_frame(4), MolError);
    ASSERT_TRUE(traj_sel.frame());
    EXPECT_EQ(traj_sel.frame(), 3);
    EXPECT_THAT(traj_sel.coords().reshaped(), ElementsAre(24, -24, 0, 48, -48, 24));
    traj_sel.coords().array() += 3;
    EXPECT_THAT(traj_sel.coords().reshaped(), ElementsAre(27, -21, 3, 51, -45, 27));
    // Undo changes
    traj_sel.coords().array() -= 3;

    // Correctly forwarding frames
    EXPECT_TRUE(traj_sel[0].frame());
    EXPECT_EQ(traj_sel[0].frame(), 3);

    /*
     * Iterators
     */
    std::vector<Atom> atoms(some_sel.begin(), some_sel.end());
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::name), {"NA", "O22", "ND"}));
    auto it = some_sel.begin();
    auto end = some_sel.end();
    EXPECT_EQ(end - it, 3);
    EXPECT_EQ((*it).name(), "NA");
    EXPECT_EQ((*++it).name(), "O22");
    it++;
    EXPECT_EQ((*it++).name(), "ND");
    EXPECT_TRUE(it == end);
    EXPECT_FALSE(it != end);

    // Correctly forwarding frames
    for (Atom atom : some_sel)
    {
        EXPECT_FALSE(atom.frame());
    }
    for (Atom atom : traj_sel)
    {
        EXPECT_TRUE(atom.frame());
    }

    /*
     * Bonds
     */
    EXPECT_THAT(all_sel.bonded()->indices(), ElementsAre(0, 2, 3));
    EXPECT_EQ(all_sel.bonded()->frame(), all_sel.frame());
    EXPECT_THAT(some_sel.bonded()->indices(), ElementsAre(0, 3));
    EXPECT_EQ(some_sel.bonded()->frame(), some_sel.frame());

    auto bonds = some_sel.bonds();
    EXPECT_EQ(bonds.size(), 1);
    EXPECT_EQ(bonds[0]->atom1(), 0);
    EXPECT_EQ(bonds[0]->atom2(), 3);

    bonds = all_sel.bonds();
    std::vector<size_t> bond_indices;
    bond_indices.reserve(4);
    for (auto b : bonds)
    {
        bond_indices.push_back(b->atom1());
        bond_indices.push_back(b->atom2());
    }
    EXPECT_THAT(bond_indices, UnorderedElementsAre(0, 2, 0, 3));

    /*
     * Conversions
     */
    EXPECT_THAT(all_sel.atoms()->indices(), ElementsAre(0, 1, 2, 3, 4, 5));
    EXPECT_EQ(all_sel.atoms()->frame(), all_sel.frame());
}

TEST(Selections, AtomSel) {
    PDBFiles pdb;
    pdb.check();

    /*
     * Construction
     */
    AtomSel all_sel(pdb.tiny);
    EXPECT_EQ(all_sel.size(), pdb.tiny->size());
    EXPECT_FALSE(all_sel.frame());
    EXPECT_THAT(all_sel.indices(), ElementsAre(0, 1, 2, 3, 4, 5));
    EXPECT_THAT(all_sel.atom_indices(), ElementsAre(0, 1, 2, 3, 4, 5));
    EXPECT_THAT(AtomSel::from_atom_indices({0, 4, 1}, pdb.tiny).indices(), ElementsAre(0, 1, 4));
    EXPECT_THAT(all_sel.selected(), ElementsAre(true, true, true, true, true, true));
    for (size_t i = 0; i < 6; ++i)
    {
        EXPECT_TRUE(all_sel.contains(i)) << "index " << i;
    }
    EXPECT_FALSE(all_sel.contains(6));

    // Constructors accepting indexes
    std::vector<size_t> indices{4, 1, 3};
    AtomSel some_sel(indices, pdb.tiny);
    AtomSel rvalue_sel({4, 1, 3}, pdb.tiny);
    EXPECT_EQ(some_sel.size(), indices.size());
    EXPECT_EQ(rvalue_sel.size(), indices.size());
    EXPECT_FALSE(some_sel.frame());
    EXPECT_FALSE(rvalue_sel.frame());
    EXPECT_THAT(some_sel.indices(), ElementsAre(1, 3, 4));
    EXPECT_THAT(rvalue_sel.indices(), ElementsAre(1, 3, 4));
    EXPECT_THAT(some_sel.atom_indices(), ElementsAre(1, 3, 4));
    EXPECT_THAT(rvalue_sel.atom_indices(), ElementsAre(1, 3, 4));
    EXPECT_THAT(some_sel.selected(), ElementsAre(false, true, false, true, true, false));
    EXPECT_THAT(rvalue_sel.selected(), ElementsAre(false, true, false, true, true, false));
    for (size_t i : {1, 3, 4})
    {
        EXPECT_TRUE(some_sel.contains(i)) << "Index " << i;
        EXPECT_TRUE(rvalue_sel.contains(i)) << "Index " << i;
    }
    for (size_t i : {0, 2, 5})
    {
        EXPECT_FALSE(some_sel.contains(i)) << "Index " << i;
        EXPECT_FALSE(rvalue_sel.contains(i)) << "Index " << i;
    }

    /*
     * Atoms getters
     */
    EXPECT_EQ(some_sel[0].index(), 1);
    EXPECT_FALSE(some_sel[0].frame());
    EXPECT_EQ(some_sel.at(0).index(), 1);
    EXPECT_FALSE(some_sel.at(0).frame());
    EXPECT_EQ(rvalue_sel[0].index(), 1);
    EXPECT_EQ(some_sel[1].index(), 3);
    EXPECT_FALSE(some_sel[1].frame());
    EXPECT_EQ(some_sel.at(1).index(), 3);
    EXPECT_FALSE(some_sel.at(1).frame());
    EXPECT_EQ(rvalue_sel[1].index(), 3);
    EXPECT_EQ(some_sel[2].index(), 4);
    EXPECT_FALSE(some_sel[2].frame());
    EXPECT_EQ(some_sel.at(2).index(), 4);
    EXPECT_FALSE(some_sel.at(2).frame());
    EXPECT_EQ(rvalue_sel[2].index(), 4);
    EXPECT_THROW(some_sel.at(3).index(), MolError);
    for (size_t i : {1, 3, 4})
    {
        EXPECT_TRUE(some_sel.contains(i)) << "Index " << i;
        EXPECT_TRUE(rvalue_sel.contains(i)) << "Index " << i;
    }
    for (size_t i : {0, 2, 5})
    {
        EXPECT_FALSE(some_sel.contains(i)) << "Index " << i;
        EXPECT_FALSE(rvalue_sel.contains(i)) << "Index " << i;
    }

    /*
     * Trajectory
     */
    EXPECT_THROW(all_sel.coords(), MolError);
    AtomSel traj_sel(pdb.traj);
    ASSERT_TRUE(traj_sel.frame());
    EXPECT_EQ(traj_sel.frame(), 0);
    EXPECT_THAT(traj_sel.coords().reshaped(), ElementsAre(1, -1, 0, 2, -2, 1));
    traj_sel.set_frame(3);
    ASSERT_TRUE(traj_sel.frame());
    EXPECT_EQ(traj_sel.frame(), 3);
    EXPECT_THROW(traj_sel.set_frame(4), MolError);
    ASSERT_TRUE(traj_sel.frame());
    EXPECT_EQ(traj_sel.frame(), 3);
    EXPECT_THAT(traj_sel.coords().reshaped(), ElementsAre(24, -24, 0, 48, -48, 24));
    traj_sel.coords().array() += 3;
    EXPECT_THAT(traj_sel.coords().reshaped(), ElementsAre(27, -21, 3, 51, -45, 27));
    // Undo changes
    traj_sel.coords().array() -= 3;

    // Correctly forwarding frames
    EXPECT_TRUE(traj_sel[0].frame());
    EXPECT_EQ(traj_sel[0].frame(), 3);

    /*
     * Iterators
     */
    std::vector<Atom> atoms(some_sel.begin(), some_sel.end());
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::name), {"NA", "O22", "ND"}));
    auto it = some_sel.begin();
    auto end = some_sel.end();
    EXPECT_EQ(end - it, 3);
    EXPECT_EQ((*it).name(), "NA");
    EXPECT_EQ((*++it).name(), "O22");
    it++;
    EXPECT_EQ((*it++).name(), "ND");
    EXPECT_TRUE(it == end);
    EXPECT_FALSE(it != end);

    // Correctly forwarding frames
    for (Atom atom : some_sel)
    {
        EXPECT_FALSE(atom.frame());
    }
    for (Atom atom : traj_sel)
    {
        EXPECT_TRUE(atom.frame());
    }

    /*
     * Bonds
     */
    EXPECT_THAT(all_sel.bonded()->indices(), ElementsAre(0, 2, 3));
    EXPECT_EQ(all_sel.bonded()->frame(), all_sel.frame());
    EXPECT_THAT(some_sel.bonded()->indices(), ElementsAre(0, 3));
    EXPECT_EQ(some_sel.bonded()->frame(), some_sel.frame());

    auto bonds = some_sel.bonds();
    EXPECT_EQ(bonds.size(), 1);
    EXPECT_EQ(bonds[0]->atom1(), 0);
    EXPECT_EQ(bonds[0]->atom2(), 3);

    bonds = all_sel.bonds();
    std::vector<size_t> bond_indices;
    bond_indices.reserve(4);
    for (auto b : bonds)
    {
        bond_indices.push_back(b->atom1());
        bond_indices.push_back(b->atom2());
    }
    EXPECT_THAT(bond_indices, UnorderedElementsAre(0, 2, 0, 3));

    /*
     * Conversions
     */
    EXPECT_THAT(all_sel.atoms()->indices(), ElementsAre(0, 1, 2, 3, 4, 5));
    EXPECT_EQ(all_sel.atoms()->frame(), all_sel.frame());
}

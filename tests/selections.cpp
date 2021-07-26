#include "files.hpp"
#include "matchers.hpp"
#include "Atom.hpp"
#include "AtomSel.hpp"
#include "MolError.hpp"
#include "internal/BaseSel.hpp"
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
    // Default constructor
    BaseSel<Atom> all_sel(pdb.tiny->size(), pdb.tiny);
    EXPECT_EQ(all_sel.size(), 6);
    EXPECT_EQ(all_sel.frame(), 0);
    EXPECT_THAT(all_sel.indices(), ElementsAre(0, 1, 2, 3, 4, 5));
    for (size_t i = 0; i < 6; ++i)
    {
        EXPECT_TRUE(all_sel.contains(i)) << "index " << i;
    }
    EXPECT_FALSE(all_sel.contains(6));

    // Constructors accepting indexes
    std::vector<size_t> indices{4, 1, 3};
    BaseSel<Atom> some_sel(pdb.tiny->size(), indices, pdb.tiny);
    BaseSel<Atom> rvalue_sel(pdb.tiny->size(), {4, 1, 3}, pdb.tiny);
    EXPECT_EQ(some_sel.size(), 3);
    EXPECT_EQ(rvalue_sel.size(), 3);
    EXPECT_EQ(some_sel.frame(), 0);
    EXPECT_EQ(rvalue_sel.frame(), 0);
    EXPECT_THAT(some_sel.indices(), ElementsAre(1, 3, 4));
    EXPECT_THAT(rvalue_sel.indices(), ElementsAre(1, 3, 4));

    /*
     * Queries
     */
    EXPECT_EQ(some_sel[0].index(), 1);
    EXPECT_EQ(some_sel.at(0).index(), 1);
    EXPECT_EQ(rvalue_sel[0].index(), 1);
    EXPECT_EQ(some_sel[1].index(), 3);
    EXPECT_EQ(some_sel.at(1).index(), 3);
    EXPECT_EQ(rvalue_sel[1].index(), 3);
    EXPECT_EQ(some_sel[2].index(), 4);
    EXPECT_EQ(some_sel.at(2).index(), 4);
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
    EXPECT_TRUE(check_valid_frame(0, pdb.traj));
    EXPECT_FALSE(check_valid_frame(4, pdb.traj));
    AtomSel traj_sel(pdb.traj->size(), pdb.traj);
    EXPECT_EQ(traj_sel[0].frame(), 0);
    traj_sel.set_frame(2);
    EXPECT_EQ(traj_sel[0].frame(), 2);
    EXPECT_THROW(traj_sel.set_frame(4), MolError);
    EXPECT_EQ(traj_sel[0].frame(), 2);

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
}

TEST(Selections, AtomSel) {
    PDBFiles pdb;
    pdb.check();

    // Default constructor
    AtomSel all_sel(pdb.tiny->size(), pdb.tiny);
    EXPECT_EQ(all_sel.size(), 6);
    AtomSel some_sel(pdb.tiny->size(), {4, 1, 3}, pdb.tiny);
    EXPECT_EQ(some_sel.size(), 3);

    /*
     * Bonds
     */
    EXPECT_THAT(all_sel.bonded()->indices(), ElementsAre(0, 2, 3));
    EXPECT_THAT(some_sel.bonded()->indices(), ElementsAre(0, 3));

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
     * Trajectory
     */
    AtomSel traj_sel(pdb.traj->size(), pdb.traj);
    EXPECT_EQ(traj_sel[0].frame(), 0);
    EXPECT_THAT(traj_sel[0].coords().reshaped(), ElementsAre(1, -1, 0));
    traj_sel.set_frame(2);
    EXPECT_THAT(traj_sel[0].coords().reshaped(), ElementsAre(6, -6, 0));
    traj_sel.set_frame(3);
    EXPECT_THAT(traj_sel.coords().reshaped(), ElementsAre(24, -24, 0, 48, -48, 24));
    traj_sel.coords().array() += 3;
    EXPECT_THAT(traj_sel.coords().reshaped(), ElementsAre(27, -21, 3, 51, -45, 27));
    // Undo changes
    traj_sel.coords().array() -= 3;
}

#include "files.hpp"
#include "matchers.hpp"
#include "Atom.hpp"
#include "AtomSel.hpp"
#include "MolError.hpp"
#include "core/AtomData.hpp"
#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace mol::internal;
using namespace testing;

TEST(Atoms, Atom) {
    auto data = AtomData::create(3);
    ASSERT_THAT(data, NotNull());

    // Comparison
    EXPECT_TRUE(Atom(1, 0, data) == Atom(1, 0, data));
    EXPECT_FALSE(Atom(0, 0, data) == Atom(1, 0, data));

    // Sizes
    ASSERT_EQ(data->size(), 3);

    /*
     * Properties
     */
    Atom atom0(0, 0, data);
    Atom atom(1, 0, data);
    EXPECT_EQ(atom.index(), 1);
    EXPECT_EQ(atom.frame(), 0);

    atom.set_resid(20);
    EXPECT_EQ(atom.resid(), 20);

    atom.set_atomic(2);
    EXPECT_EQ(atom.atomic(), 2);

    atom.set_occupancy(1.0);
    EXPECT_EQ(atom.occupancy(), 1.0);

    atom.set_tempfactor(1.0);
    EXPECT_EQ(atom.tempfactor(), 1.0);

    atom.set_mass(14.0);
    EXPECT_EQ(atom.mass(), 14.0);

    atom.set_charge(-1.0);
    EXPECT_EQ(atom.charge(), -1.0);

    atom.set_radius(1.0);
    EXPECT_EQ(atom.radius(), 1.0);

    atom.set_name("CA");
    EXPECT_EQ(atom.name(), "CA");

    atom.set_type("C");
    EXPECT_EQ(atom.type(), "C");

    atom.set_resname("ARG");
    EXPECT_EQ(atom.resname(), "ARG");

    atom.set_segid("SEG1");
    EXPECT_EQ(atom.segid(), "SEG1");

    atom.set_chain("B");
    EXPECT_EQ(atom.chain(), "B");

    atom.set_altloc("B");
    EXPECT_EQ(atom.altloc(), "B");

    /*
     * Coordinates
     */
    Timestep ts(3);
    ts.coords() << 1, 2, 3, 4, 5, 6, 7, 8, 9;

    ASSERT_EQ(data->num_frames(), 0);
    data->add_timestep(std::move(ts));
    ASSERT_EQ(data->num_frames(), 1);

    EXPECT_THAT(atom.coords(), ElementsAre(2, 5, 8));
    atom.coords() *= 2;
    EXPECT_THAT(atom.coords(), ElementsAre(4, 10, 16));

    /*
     * Bonds
     */
    ASSERT_THAT(atom.add_bond(2), NotNull());
    EXPECT_THROW(atom.add_bond(1), MolError);
    EXPECT_THROW(atom.add_bond(3), MolError);

    EXPECT_THAT(atom0.bond(2), IsNull());
    EXPECT_THAT(atom.bond(3), IsNull());
    auto bond = atom.bond(2);
    ASSERT_THAT(bond, NotNull());
    EXPECT_EQ(bond->atom1(), 1);
    EXPECT_EQ(bond->atom2(), 2);

    EXPECT_EQ(atom0.bonds().size(), 0);
    auto bond_list = atom.bonds();
    EXPECT_EQ(bond_list.size(), 1);
    EXPECT_EQ(bond_list[0]->atom1(), 1);
    EXPECT_EQ(bond_list[0]->atom2(), 2);

    EXPECT_EQ(atom0.bonded()->size(), 0);
    auto bond_sel = atom.bonded();
    EXPECT_EQ(bond_sel->size(), 1);
    EXPECT_THAT(bond_sel->indices(), ElementsAre(2));
}

TEST(Atoms, Timestep) {
    Timestep ts(2);

    ts.coords() << 1, 2, 3, 4, 5, 6;
    EXPECT_THAT(ts.coords().reshaped(), ElementsAre(1, 3, 5, 2, 4, 6));

    // Note: ts is invalid from this on
    float *data = ts.coords().data();
    Timestep moved(std::move(ts));
    EXPECT_THAT(ts.coords().data(), IsNull());
    ASSERT_EQ(moved.coords().data(), data);
    EXPECT_THAT(moved.coords().reshaped(), ElementsAre(1, 3, 5, 2, 4, 6));

    Timestep moved_again = std::move(moved);
    EXPECT_THAT(moved.coords().data(), IsNull());
    ASSERT_EQ(moved_again.coords().data(), data);
    EXPECT_THAT(moved_again.coords().reshaped(), ElementsAre(1, 3, 5, 2, 4, 6));
}

TEST(Atoms, AtomSel) {
    PDBFiles pdb;
    pdb.check();

    // Default constructor
    AtomSel all_sel(pdb.tiny);
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
    AtomSel some_sel(indices, pdb.tiny);
    AtomSel rvalue_sel({4, 1, 3}, pdb.tiny);
    EXPECT_EQ(some_sel.size(), 3);
    EXPECT_EQ(rvalue_sel.size(), 3);
    EXPECT_EQ(some_sel.frame(), 0);
    EXPECT_EQ(rvalue_sel.frame(), 0);
    EXPECT_THAT(some_sel.indices(), ElementsAre(1, 3, 4));
    EXPECT_THAT(rvalue_sel.indices(), ElementsAre(1, 3, 4));
    EXPECT_EQ(some_sel[0].index(), 1);
    EXPECT_EQ(rvalue_sel[0].index(), 1);
    EXPECT_EQ(some_sel[1].index(), 3);
    EXPECT_EQ(rvalue_sel[1].index(), 3);
    EXPECT_EQ(some_sel[2].index(), 4);
    EXPECT_EQ(rvalue_sel[2].index(), 4);
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
    AtomSel traj_sel(pdb.traj);
    EXPECT_EQ(traj_sel[0].frame(), 0);
    EXPECT_THAT(traj_sel[0].coords().reshaped(), ElementsAre(1, -1, 0));
    traj_sel.set_frame(2);
    EXPECT_EQ(traj_sel[0].frame(), 2);
    EXPECT_THAT(traj_sel[0].coords().reshaped(), ElementsAre(6, -6, 0));
    EXPECT_THROW(traj_sel.set_frame(4), MolError);
    EXPECT_EQ(traj_sel[0].frame(), 2);

    /*
     * Iterators
     */
    std::vector<Atom> atoms(some_sel.begin(), some_sel.end());
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::resid), {339, 801, 85}));
    auto it = some_sel.begin();
    auto end = some_sel.end();
    EXPECT_EQ(end - it, 3);
    EXPECT_EQ((*it).resid(), 339);
    EXPECT_EQ((*++it).resid(), 801);
    it++;
    EXPECT_EQ((*it++).resid(), 85);
    EXPECT_TRUE(it == end);
    EXPECT_FALSE(it != end);

    /*
     * Accessors
     */
    traj_sel.set_frame(3);
    EXPECT_THAT(traj_sel.coords().reshaped(), ElementsAre(24, -24, 0, 48, -48, 24));
    traj_sel.coords().array() += 3;
    EXPECT_THAT(traj_sel.coords().reshaped(), ElementsAre(27, -21, 3, 51, -45, 27));
}

TEST(Atoms, bonds) {
    PDBFiles pdb;
    pdb.check();
    AtomSel all_sel(pdb.tiny);

    EXPECT_EQ(all_sel[1].bonded()->size(), 0);
    EXPECT_EQ(all_sel[4].bonded()->size(), 0);
    EXPECT_EQ(all_sel[5].bonded()->size(), 0);
    EXPECT_THAT(all_sel[0].bonded()->indices(), ElementsAre(2, 3));
    EXPECT_THAT(all_sel[2].bonded()->indices(), ElementsAre(0));
    EXPECT_THAT(all_sel[3].bonded()->indices(), ElementsAre(0));
}

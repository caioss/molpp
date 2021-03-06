#include "files.hpp"
#include "matchers.hpp"
#include "auxiliary.hpp"
#include "core/MolData.hpp"
#include "core/AtomData.hpp"
#include <molpp/Atom.hpp>
#include <molpp/Residue.hpp>
#include <molpp/AtomSel.hpp>
#include <molpp/MolError.hpp>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace mol::internal;
using namespace testing;

TEST(Atoms, Atom) {
    auto data = create_moldata(3, 1, 1, 1, 1);

    // Comparison
    EXPECT_TRUE(Atom(1, 0, data) == Atom(1, 0, data));
    EXPECT_FALSE(Atom(0, 0, data) == Atom(1, 0, data));
    EXPECT_FALSE(Atom(1, std::nullopt, data) == Atom(1, 0, data));
    EXPECT_FALSE(Atom(1, 0, data) == Atom(1, 0, nullptr));

    /*
     * Properties
     */
    Atom atom0(0, std::nullopt, data);
    Atom atom(1, 0, data);
    EXPECT_EQ(atom.index(), 1);
    EXPECT_EQ(atom.frame(), 0);
    EXPECT_EQ(atom0.frame(), std::nullopt);

    EXPECT_EQ(atom.resid(), 1);

    EXPECT_EQ(atom.residue_id(), 1);
    EXPECT_EQ(atom.residue(), Residue(1, 0, data));

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

    EXPECT_EQ(atom.resname(), "B");
    EXPECT_EQ(atom.segid(), "A");
    EXPECT_EQ(atom.chain(), "A");

    atom.set_altloc("B");
    EXPECT_EQ(atom.altloc(), "B");

    /*
     * Coordinates
     */
    EXPECT_THAT(atom.coords().reshaped(), ElementsAre(1, 1, 1));
    atom.coords() *= 2;
    EXPECT_THAT(atom.coords().reshaped(), ElementsAre(2, 2, 2));

    EXPECT_THROW(Atom(1, std::nullopt, data).coords(), MolError);

    /*
     * Bonds
     */
    ASSERT_THAT(atom.add_bond(2), NotNull());
    EXPECT_THROW(atom.add_bond(1), MolError);
    EXPECT_THROW(atom.add_bond(3), MolError);
    EXPECT_EQ(atom.add_bond(2), atom.bond(2));
    EXPECT_EQ(atom.add_bond(Atom(2, 0, data)), atom.bond(Atom(2, 0, data)));

    EXPECT_THAT(atom0.bond(1), NotNull());
    EXPECT_THAT(atom.bond(2), NotNull());
    auto bond = atom.bond(2);
    ASSERT_THAT(bond, NotNull());
    EXPECT_EQ(bond->atom1(), 1);
    EXPECT_EQ(bond->atom2(), 2);

    EXPECT_EQ(atom0.bonds().size(), 1);
    auto bond_list = atom0.bonds();
    ASSERT_EQ(bond_list.size(), 1);
    EXPECT_EQ(bond_list[0]->atom1(), 0);
    EXPECT_EQ(bond_list[0]->atom2(), 1);
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

TEST(Atoms, MolData) {
    size_t const num_atoms { 3 };
    auto data = MolData::create(num_atoms);
    ASSERT_THAT(data, NotNull());
    EXPECT_EQ(data->size(), num_atoms);
    EXPECT_EQ(data->atoms().size(), num_atoms);
    EXPECT_EQ(data->bonds().size(), 0);
    EXPECT_EQ(data->residues().size(), 0);
    EXPECT_EQ(data->trajectory().num_frames(), 0);
}

TEST(Atoms, Trajectory) {
    size_t const num_atoms { 3 };
    auto data = MolData::create(num_atoms);
    ASSERT_THAT(data, NotNull());

    Trajectory &traj_data = data->trajectory();
    EXPECT_EQ(traj_data.num_frames(), 0);
    traj_data.add_timestep(Timestep(num_atoms));
    EXPECT_EQ(traj_data.num_frames(), 1);
    EXPECT_EQ(traj_data.timestep(0).coords().cols(), num_atoms);
}

TEST(Atoms, AtomData) {
    AtomData props(1);
    EXPECT_EQ(props.size(), 1);
    EXPECT_EQ(props.residue(0), -1);
    EXPECT_EQ(props.atomic(0), 0);
    EXPECT_EQ(props.occupancy(0), 0);
    EXPECT_EQ(props.tempfactor(0), 0);
    EXPECT_EQ(props.mass(0), 0);
    EXPECT_EQ(props.charge(0), 0);
    EXPECT_EQ(props.radius(0), 0);
    EXPECT_EQ(props.name(0), "");
    EXPECT_EQ(props.type(0), "");
    EXPECT_EQ(props.altloc(0), "");
}

#include "files.hpp"
#include "matchers.hpp"
#include "auxiliary.hpp"
#include <molpp/internal/AtomData.hpp>
#include <molpp/internal/MolData.hpp>
#include <molpp/Atom.hpp>
#include <molpp/Residue.hpp>
#include <molpp/AtomSel.hpp>
#include <molpp/MolError.hpp>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace mol::internal;
using namespace testing;

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
    MolData data(num_atoms);
    EXPECT_EQ(data.size(), num_atoms);
    EXPECT_EQ(data.atoms().size(), num_atoms);
    EXPECT_EQ(data.bonds().size(), 0);
    EXPECT_EQ(data.residues().size(), 0);
    EXPECT_EQ(data.trajectory().num_frames(), 0);
}

TEST(Atoms, Trajectory) {
    size_t const num_atoms { 3 };
    MolData data(num_atoms);
    EXPECT_EQ(data.size(), num_atoms);

    Trajectory &traj_data = data.trajectory();
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

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "Atom.hpp"
#include "core/AtomData.hpp"

using namespace testing;
using namespace mol::internal;

TEST(atomdata, BasicAssertions) {
    auto data = AtomData::create(2);
    auto atom = data->index(1);

    // Sizes
    ASSERT_EQ(data->size(), 2);

    // Properties
    ASSERT_THAT(data->indexes(), ElementsAre(0, 1));
    EXPECT_EQ(atom.index(), 1);

    atom.set_resid(20);
    EXPECT_EQ(atom.resid(), 20);
    ASSERT_THAT(data->resids(), ElementsAre(-1, 20));

    atom.set_atomic(2);
    EXPECT_EQ(atom.atomic(), 2);
    ASSERT_THAT(data->atomics(), ElementsAre(0, 2));

    atom.set_occupancy(1.0);
    EXPECT_EQ(atom.occupancy(), 1.0);
    ASSERT_THAT(data->occupancies(), ElementsAre(0, 1.0));

    atom.set_tempfactor(1.0);
    EXPECT_EQ(atom.tempfactor(), 1.0);
    ASSERT_THAT(data->tempfactors(), ElementsAre(0, 1.0));

    atom.set_mass(14.0);
    EXPECT_EQ(atom.mass(), 14.0);
    ASSERT_THAT(data->masses(), ElementsAre(0, 14.0));

    atom.set_charge(-1.0);
    EXPECT_EQ(atom.charge(), -1.0);
    ASSERT_THAT(data->charges(), ElementsAre(0, -1.0));

    atom.set_radius(1.0);
    EXPECT_EQ(atom.radius(), 1.0);
    ASSERT_THAT(data->radii(), ElementsAre(0, 1.0));

    atom.set_name("CA");
    EXPECT_EQ(atom.name(), "CA");
    ASSERT_THAT(data->names(), ElementsAre("", "CA"));

    atom.set_type("C");
    EXPECT_EQ(atom.type(), "C");
    ASSERT_THAT(data->types(), ElementsAre("", "C"));

    atom.set_resname("ARG");
    EXPECT_EQ(atom.resname(), "ARG");
    ASSERT_THAT(data->resnames(), ElementsAre("", "ARG"));

    atom.set_segid("SEG1");
    EXPECT_EQ(atom.segid(), "SEG1");
    ASSERT_THAT(data->segids(), ElementsAre("", "SEG1"));

    atom.set_chain("B");
    EXPECT_EQ(atom.chain(), "B");
    ASSERT_THAT(data->chains(), ElementsAre("", "B"));

    atom.set_altloc("B");
    EXPECT_EQ(atom.altloc(), "B");
    ASSERT_THAT(data->altlocs(), ElementsAre("", "B"));

    // TODO
    // TODO add coords accesssors to atoms
}

TEST(timestep, BasicAssertions) {
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

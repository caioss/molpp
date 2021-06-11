#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "Atom.hpp"
#include "core/AtomData.hpp"

using namespace testing;
using namespace mol::internal;

TEST(atoms, BasicAssertions) {
    auto data = AtomData::create(2);
    ASSERT_THAT(data, NotNull());
    auto atom = data->index(1, 0);

    // Comparison
    EXPECT_TRUE(data->index(1, 0) == data->index(1, 0));
    EXPECT_FALSE(data->index(0, 0) == data->index(1, 0));

    // Sizes
    ASSERT_EQ(data->size(), 2);

    /*
     * Properties
     */
    EXPECT_EQ(atom.index(), 1);

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
    Timestep ts(2);
    ts.coords() << 1, 2, 3, 4, 5, 6;

    ASSERT_EQ(data->num_frames(), 0);
    data->add_timestep(std::move(ts));
    ASSERT_EQ(data->num_frames(), 1);

    EXPECT_THAT(atom.coords(), ElementsAre(2, 4, 6));
    atom.coords() *= 2;
    EXPECT_THAT(atom.coords(), ElementsAre(4, 8, 12));
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

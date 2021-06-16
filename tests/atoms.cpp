#include "matchers.hpp"
#include "Atom.hpp"
#include "AtomSel.hpp"
#include "MolError.hpp"
#include "core/AtomData.hpp"
#include "readers/MolReader.hpp"
#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace mol::internal;
using namespace testing;

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

TEST(atomsel, BasicAssertions) {
    auto reader = MolReader::from_file_ext(".pdb");
    auto atom_data = reader->read_topology("tiny.pdb");
    reader->read_trajectory("tiny.pdb", atom_data);

    AtomSel all_sel(atom_data);
    EXPECT_EQ(all_sel.size(), 6);
    EXPECT_EQ(all_sel.frame(), 0);
    EXPECT_THAT(all_sel.indices(), ElementsAre(0, 1, 2, 3, 4, 5));
    for (size_t i = 0; i < 6; ++i)
    {
        EXPECT_TRUE(all_sel.contains(i)) << "index " << i;
    }
    EXPECT_FALSE(all_sel.contains(6));

    AtomSel messy_sel({3, 1, 4, 3}, atom_data);
    EXPECT_EQ(messy_sel.size(), 3);
    EXPECT_EQ(messy_sel.frame(), 0);
    EXPECT_THAT(messy_sel.indices(), ElementsAre(1, 3, 4));
    for (size_t i : {1, 3, 4})
    {
        EXPECT_TRUE(messy_sel.contains(i)) << "Index " << i;
    }
    for (size_t i : {0, 2, 5})
    {
        EXPECT_FALSE(messy_sel.contains(i)) << "Index " << i;
    }

    /*
     * Trajectory
     */
    auto traj_data = reader->read_topology("traj.pdb");
    reader->read_trajectory("traj.pdb", traj_data);
    AtomSel traj_sel(traj_data);
    EXPECT_THAT(traj_sel[0].coords().reshaped(), ElementsAre(1, -1, 0));
    traj_sel.set_frame(2);
    EXPECT_THAT(traj_sel[0].coords().reshaped(), ElementsAre(6, -6, 0));
    EXPECT_THROW(traj_sel.set_frame(4), MolError);

    /*
     * Iterators
     */
    std::vector<Atom> atoms(messy_sel.begin(), messy_sel.end());
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::resid), {339, 801, 85}));
    auto it = messy_sel.begin();
    auto end = messy_sel.end();
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

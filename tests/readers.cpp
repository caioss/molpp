#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "MolError.hpp"
#include "core/AtomData.hpp"
#include "readers/MolReader.hpp"
#include "readers/MolfileReader.hpp"

using namespace testing;
using namespace mol::internal;

TEST(molfile_structure, BasicAssertions) {
    EXPECT_THROW(MolfileReader(""), mol::MolError);

    MolfileReader reader("pdb");
    ASSERT_TRUE(reader.open("tiny.pdb"));
    auto atoms = reader.read_atoms();
    reader.close();

    // Size
    ASSERT_THAT(atoms, NotNull());
    ASSERT_EQ(atoms->size(), 6);

    // Properties
    ASSERT_THAT(atoms->resids(), ElementsAre(3, 339, 201, 801, 85, 85));
    ASSERT_THAT(atoms->atomics(), ElementsAre(7, 7, 8, 8, 7, 7));
    ASSERT_THAT(atoms->occupancies(), ElementsAre(999, 1, 1, 1, 1, 1));
    ASSERT_THAT(atoms->tempfactors(), ElementsAre(-99, -1, -1, -1, -1, -1));
    ASSERT_THAT(atoms->masses(), ElementsAre(FloatEq(14.0067), FloatEq(14.0067), FloatEq(15.9994), FloatEq(15.9994), FloatEq(14.0067), FloatEq(14.0067)));
    ASSERT_THAT(atoms->radii(), ElementsAre(FloatEq(1.55), FloatEq(1.55), FloatEq(1.52), FloatEq(1.52), FloatEq(1.55), FloatEq(1.55)));
    ASSERT_THAT(atoms->names(), ElementsAre("N", "NA", "O1", "O22", "ND", "ND"));
    ASSERT_THAT(atoms->types(), ElementsAre("N", "NA", "O1", "O22", "ND", "ND"));
    ASSERT_THAT(atoms->resnames(), ElementsAre("GLY", "ASP", "HOH", "HOH", "ASP", "ASP"));
    ASSERT_THAT(atoms->segids(), ElementsAre("SEG1", "SEG2", "SEG1", "SEG2", "SEG1", "SEG1"));
    ASSERT_THAT(atoms->chains(), ElementsAre("A", "B", "A", "B", "A", "A"));
    ASSERT_THAT(atoms->altlocs(), ElementsAre(" ", " ", " ", " ", "A", "B"));

    // TODO add test to charges
}

TEST(molfile_trajectory, BasicAssertions) {
    MolfileReader reader("pdb");
    ASSERT_TRUE(reader.open("traj.pdb"));
    auto atoms = reader.read_atoms();
    reader.close();

    // Size
    ASSERT_THAT(atoms, NotNull());
    ASSERT_EQ(atoms->size(), 2);

    // Trajectory
    ASSERT_EQ(atoms->num_frames(), 0);
    atoms->add_timestep(Timestep(atoms->size()));
    ASSERT_EQ(atoms->num_frames(), 1);

    ASSERT_TRUE(reader.open("traj.pdb"));

    ASSERT_TRUE(reader.skip_timestep(atoms));
    ASSERT_EQ(atoms->num_frames(), 1);

    ASSERT_TRUE(reader.read_timestep(atoms));
    ASSERT_EQ(atoms->num_frames(), 2);
    EXPECT_THAT(atoms->timestep(1).coords().reshaped(), ElementsAre(2, -2, 0, 4, -4, 2));

    ASSERT_TRUE(reader.skip_timestep(atoms));
    ASSERT_EQ(atoms->num_frames(), 2);

    ASSERT_TRUE(reader.read_timestep(atoms));
    ASSERT_EQ(atoms->num_frames(), 3);
    EXPECT_THAT(atoms->timestep(2).coords().reshaped(), ElementsAre(24, -24, 0, 48, -48, 24));

    ASSERT_FALSE(reader.read_timestep(atoms));
    ASSERT_FALSE(reader.skip_timestep(atoms));

    reader.close();
}

TEST(PDB_plugin_types, BasicAssertions) {
    MolfileReader reader("pdb");
    EXPECT_TRUE(reader.has_topology());
    EXPECT_TRUE(reader.has_trajectory());
    EXPECT_FALSE(reader.has_trajectory_metadata());
    EXPECT_TRUE(reader.has_bonds());
    EXPECT_TRUE(reader.can_read("4lad.pdb"));
    EXPECT_TRUE(reader.can_read("dummy.ent"));
    EXPECT_TRUE(reader.can_read("dummy.dummy.pdb"));
    EXPECT_FALSE(reader.can_read("dummy.psf"));
    EXPECT_FALSE(reader.can_read("dummy.pd"));
    EXPECT_FALSE(reader.can_read("dummypdb"));
}

TEST(molreader, BasicAssertions) {
    EXPECT_THAT(MolReader::from_file("dummy.unk"), IsNull());
    EXPECT_THAT(MolReader::from_file("dummy.pdb"), NotNull());

    auto pdb_reader = MolReader::from_file("traj.pdb");
    auto atoms = pdb_reader->read_topology("traj.pdb");
    ASSERT_THAT(atoms, NotNull());
    EXPECT_EQ(atoms->size(), 2);
    // Fast check. The complete reading test is in the plugins tests
    EXPECT_THAT(atoms->names(), ElementsAre("N", "CA"));

    EXPECT_TRUE(pdb_reader->read_trajectory("traj.pdb", atoms));
    EXPECT_EQ(atoms->num_frames(), 4);
    EXPECT_TRUE(pdb_reader->read_trajectory("traj.pdb", atoms, 0, 1, 0));
    EXPECT_EQ(atoms->num_frames(), 5);
    EXPECT_TRUE(pdb_reader->read_trajectory("traj.pdb", atoms, 1, 2, 2));
    EXPECT_EQ(atoms->num_frames(), 6);
    EXPECT_TRUE(pdb_reader->read_trajectory("traj.pdb", atoms, 2, 0));
    EXPECT_TRUE(pdb_reader->read_trajectory("traj.pdb", atoms, -2, 0, 2));
    EXPECT_EQ(atoms->num_frames(), 6);
}

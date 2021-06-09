#include "matchers.hpp"
#include "MolError.hpp"
#include "Atom.hpp"
#include "core/AtomData.hpp"
#include "readers/MolReader.hpp"
#include "readers/MolfileReader.hpp"
#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace testing;
using namespace mol;
using namespace mol::internal;

TEST(molfile_structure, BasicAssertions) {
    EXPECT_THROW(MolfileReader(""), mol::MolError);

    MolfileReader reader("pdb");
    ASSERT_TRUE(reader.open("tiny.pdb"));
    auto data = reader.read_atoms();
    reader.close();

    // Size
    ASSERT_THAT(data, NotNull());
    ASSERT_EQ(data->size(), 6);

    // Properties
    std::vector<mol::Atom> atoms;
    for (size_t i = 0; i < data->size(); ++i)
    {
        atoms.push_back(data->index(i));
    }
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::resid),
                                 {3, 339, 201, 801, 85, 85}));
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::atomic),
                                 {7, 7, 8, 8, 7, 7}));
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::occupancy),
                                 {999, 1, 1, 1, 1, 1}));
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::tempfactor),
                                 {-99, -1, -1, -1, -1, -1}));
    EXPECT_THAT(atoms, Pointwise(PropFloat(&Atom::mass, 1e-5),
                                 {14.0067, 14.0067, 15.9994, 15.9994, 14.0067, 14.0067}));
    EXPECT_THAT(atoms, Pointwise(PropFloat(&Atom::radius, 1e-5),
                                 {1.55, 1.55, 1.52, 1.52, 1.55, 1.55}));
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::name),
                                 {"N", "NA", "O1", "O22", "ND", "ND"}));
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::type),
                                 {"N", "NA", "O1", "O22", "ND", "ND"}));
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::resname),
                                 {"GLY", "ASP", "HOH", "HOH", "ASP", "ASP"}));
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::segid),
                                 {"SEG1", "SEG2", "SEG1", "SEG2", "SEG1", "SEG1"}));
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::chain),
                                 {"A", "B", "A", "B", "A", "A"}));
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::altloc),
                                 {" ", " ", " ", " ", "A", "B"}));

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
    EXPECT_EQ(atoms->index(0).name(), "N");
    EXPECT_EQ(atoms->index(1).name(), "CA");

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

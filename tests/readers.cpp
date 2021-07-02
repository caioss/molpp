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

TEST(Readers, MolfileReader) {
    EXPECT_THROW(MolfileReader(""), MolError);

    MolfileReader reader(".pdb");

    // Sanity checks
    EXPECT_THROW(reader.read_atoms(), MolError); // Not open
    ASSERT_EQ(reader.open("tiny.pdb"), MolReader::SUCCESS);
    ASSERT_EQ(reader.open("tiny.pdb"), MolReader::INVALID); // Re-open

    auto data = reader.read_atoms();
    reader.close();

    ASSERT_THAT(data, NotNull());
    ASSERT_EQ(data->size(), 6);

    /*
     * Properties
     */
    std::vector<mol::Atom> atoms;
    for (size_t i = 0; i < data->size(); ++i)
    {
        atoms.push_back(data->index(i, 0));
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

    /*
     * Trajectory
     */
    ASSERT_EQ(reader.open("traj.pdb"), MolReader::SUCCESS);
    data = reader.read_atoms();
    reader.close();
    ASSERT_THAT(data, NotNull());
    ASSERT_EQ(data->size(), 2);

    // Sanity checks
    reader.open("tiny.pdb"); // Wrong number of atoms
    EXPECT_EQ(reader.check_timestep_read(data), MolReader::WRONG_ATOMS);
    reader.close(); // No file handle
    EXPECT_EQ(reader.check_timestep_read(data), MolReader::INVALID);

    reader.open("traj.pdb");

    ASSERT_EQ(reader.skip_timestep(data), MolReader::SUCCESS);
    ASSERT_EQ(data->num_frames(), 0);

    ASSERT_EQ(reader.read_timestep(data), MolReader::SUCCESS);
    ASSERT_EQ(data->num_frames(), 1);
    EXPECT_THAT(data->timestep(0).coords().reshaped(), ElementsAre(2, -2, 0, 4, -4, 2));

    ASSERT_EQ(reader.skip_timestep(data), MolReader::SUCCESS);
    ASSERT_EQ(data->num_frames(), 1);

    ASSERT_EQ(reader.read_timestep(data), MolReader::SUCCESS);
    ASSERT_EQ(data->num_frames(), 2);
    EXPECT_THAT(data->timestep(1).coords().reshaped(), ElementsAre(24, -24, 0, 48, -48, 24));

    ASSERT_EQ(reader.read_timestep(data), MolReader::END);
    ASSERT_EQ(reader.skip_timestep(data), MolReader::END);

    reader.close();
}

TEST(Readers, PDB) {
    MolfileReader reader(".pdb");
    EXPECT_TRUE(reader.has_topology());
    EXPECT_TRUE(reader.has_trajectory());
    EXPECT_FALSE(reader.has_trajectory_metadata());
    EXPECT_TRUE(reader.has_bonds());
    EXPECT_TRUE(reader.can_read(".pdb"));
    EXPECT_TRUE(reader.can_read(".ent"));
    EXPECT_FALSE(reader.can_read(".psf"));
    EXPECT_FALSE(reader.can_read(".pd"));
    EXPECT_FALSE(reader.can_read(""));
}

TEST(Readers, MolReader) {
    EXPECT_THAT(MolReader::from_file_ext(".unk"), IsNull());
    EXPECT_THAT(MolReader::from_file_ext(".pdb"), NotNull());

    auto pdb_reader = MolReader::from_file_ext(".pdb");
    auto atoms = pdb_reader->read_topology("traj.pdb");
    ASSERT_THAT(atoms, NotNull());
    EXPECT_EQ(atoms->size(), 2);
    // Fast check. The complete reading test is in the plugins tests
    EXPECT_EQ(atoms->index(0, 0).name(), "N");
    EXPECT_EQ(atoms->index(1, 0).name(), "CA");

    // Sanity checks
    // EXPECT_EQ(MolReader::from_file_ext(".psf")->read_trajectory("traj.pdb", atoms), MolReader::INVALID); // TODO enable when PSF reader becomes available
    EXPECT_EQ(pdb_reader->read_trajectory("", atoms), MolReader::FAILED);
    EXPECT_EQ(pdb_reader->read_trajectory("tiny.pdb", atoms), MolReader::WRONG_ATOMS);

    EXPECT_EQ(pdb_reader->read_trajectory("traj.pdb", atoms), MolReader::SUCCESS);
    EXPECT_EQ(atoms->num_frames(), 4);
    EXPECT_EQ(pdb_reader->read_trajectory("traj.pdb", atoms, 0, 1, 0), MolReader::SUCCESS);
    EXPECT_EQ(atoms->num_frames(), 5);
    EXPECT_EQ(pdb_reader->read_trajectory("traj.pdb", atoms, 1, 2, 2), MolReader::SUCCESS);
    EXPECT_EQ(atoms->num_frames(), 6);
    EXPECT_EQ(pdb_reader->read_trajectory("traj.pdb", atoms, 2, 0), MolReader::SUCCESS);
    EXPECT_EQ(pdb_reader->read_trajectory("traj.pdb", atoms, -2, 0, 2), MolReader::SUCCESS);
    EXPECT_EQ(atoms->num_frames(), 6);

    auto atom = atoms->index(0, 0);
}

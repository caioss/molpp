#include "matchers.hpp"
#include "readers/MolReader.hpp"
#include "readers/MolfileReader.hpp"
#include "core/MolData.hpp"
#include <molpp/MolError.hpp>
#include <molpp/Atom.hpp>
#include <molpp/AtomSel.hpp>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace testing;
using namespace mol;
using namespace mol::internal;

TEST(Readers, MolfileReader) {
    EXPECT_THROW(MolfileReader(""), MolError);

    /*
     * Readers
     */
    MolfileReader reader(".pdb");

    // General checks
    EXPECT_THROW(reader.read_atoms(), MolError); // Not open
    ASSERT_EQ(reader.open("tiny.pdb"), MolReader::SUCCESS);
    ASSERT_EQ(reader.open("tiny.pdb"), MolReader::INVALID); // Re-open

    // Read PDB
    auto data = reader.read_atoms();
    reader.close();
    ASSERT_THAT(data, NotNull());
    ASSERT_EQ(data->size(), 6);

    // Read mol2
    reader = MolfileReader(".mol2");
    ASSERT_EQ(reader.open("fluorobenzene.mol2"), MolReader::SUCCESS);
    auto m2_data = reader.read_atoms();
    reader.close();
    ASSERT_THAT(m2_data, NotNull());
    ASSERT_EQ(m2_data->size(), 12);

    /*
     * Atom properties
     */
    std::vector<mol::Atom> atoms;
    for (size_t i = 0; i < data->size(); ++i)
    {
        atoms.push_back(Atom(i, 0, data));
    }
    std::vector<mol::Atom> m2_atoms;
    for (size_t i = 0; i < m2_data->size(); ++i)
    {
        m2_atoms.push_back(Atom(i, 0, m2_data));
    }
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::resid),
                                 {3, 339, 201, 801, 85, 85}));
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::residue_id),
                                 {0, 1, 2, 3, 4, 4}));
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
    EXPECT_THAT(m2_atoms, Pointwise(PropFloat(&Atom::charge, 1e-5),
                                    {-0.2055,  0.1234, -0.0265, -0.0265,
                                     -0.0590, -0.0590, -0.0616,  0.0646,
                                      0.0646,  0.0618,  0.0618,  0.0618}));

    /*
     * Trajectory
     */
    reader = MolfileReader(".pdb");
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

    /*
     * Bonds and charges
     */
    // Bond partners
    EXPECT_THAT(m2_atoms[0].bonded()->indices(), ElementsAre(1));
    EXPECT_THAT(m2_atoms[1].bonded()->indices(), ElementsAre(0, 2, 3));
    EXPECT_THAT(m2_atoms[2].bonded()->indices(), ElementsAre(1, 4, 7));
    EXPECT_THAT(m2_atoms[3].bonded()->indices(), ElementsAre(1, 5, 8));
    EXPECT_THAT(m2_atoms[4].bonded()->indices(), ElementsAre(2, 6, 9));
    EXPECT_THAT(m2_atoms[5].bonded()->indices(), ElementsAre(3, 6, 10));
    EXPECT_THAT(m2_atoms[6].bonded()->indices(), ElementsAre(4, 5, 11));
    EXPECT_THAT(m2_atoms[7].bonded()->indices(), ElementsAre(2));
    EXPECT_THAT(m2_atoms[8].bonded()->indices(), ElementsAre(3));
    EXPECT_THAT(m2_atoms[9].bonded()->indices(), ElementsAre(4));
    EXPECT_THAT(m2_atoms[10].bonded()->indices(), ElementsAre(5));
    EXPECT_THAT(m2_atoms[11].bonded()->indices(), ElementsAre(6));

    // Bond orders
    for (Atom &atom : m2_atoms)
    {
        for (auto bond : atom.bonds())
        {
            EXPECT_FALSE(bond->guessed());
            EXPECT_FALSE(bond->guessed_order());
        }
    }
    EXPECT_EQ(m2_atoms[0].bond(1)->order(), Bond::Single);
    EXPECT_EQ(m2_atoms[1].bond(2)->order(), Bond::Aromatic);
    EXPECT_EQ(m2_atoms[1].bond(3)->order(), Bond::Aromatic);
    EXPECT_EQ(m2_atoms[2].bond(4)->order(), Bond::Aromatic);
    EXPECT_EQ(m2_atoms[2].bond(7)->order(), Bond::Single);
    EXPECT_EQ(m2_atoms[3].bond(5)->order(), Bond::Aromatic);
    EXPECT_EQ(m2_atoms[3].bond(8)->order(), Bond::Single);
    EXPECT_EQ(m2_atoms[4].bond(6)->order(), Bond::Aromatic);
    EXPECT_EQ(m2_atoms[4].bond(9)->order(), Bond::Single);
    EXPECT_EQ(m2_atoms[5].bond(6)->order(), Bond::Aromatic);
    EXPECT_EQ(m2_atoms[5].bond(10)->order(), Bond::Single);
    EXPECT_EQ(m2_atoms[6].bond(11)->order(), Bond::Single);

    /*
     * Residue detection
     */
    reader = MolfileReader(".pdb");
    ASSERT_EQ(reader.open("4lad.pdb"), MolReader::SUCCESS);
    data = reader.read_atoms();
    reader.close();
    ASSERT_THAT(data, NotNull());
    EXPECT_EQ(data->properties().residue(0), 0);
    EXPECT_EQ(data->properties().residue(3), 0);
    EXPECT_EQ(data->properties().residue(4), 1);
    EXPECT_EQ(data->properties().residue(1230), 153);
    EXPECT_EQ(data->properties().residue(1231), 154);
    EXPECT_EQ(data->properties().residue(1561), 195);
    EXPECT_EQ(data->properties().residue(1562), 196);
    EXPECT_EQ(data->properties().residue(1790), 222);
    EXPECT_EQ(data->properties().residue(1791), 223);
    EXPECT_EQ(data->properties().residue(1792), 224);
    EXPECT_EQ(data->properties().residue(1793), 225);
    EXPECT_EQ(data->properties().residue(1798), 225);
    EXPECT_EQ(data->properties().residue(1799), 226);
    EXPECT_EQ(data->properties().residue(1804), 226);
    // Water molecules
    for (size_t i = 0; i < 39; ++i)
    {
        EXPECT_EQ(data->properties().residue(1805 + i), 227 + i);
    }
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

TEST(Readers, Mol2) {
    MolfileReader reader(".mol2");
    EXPECT_TRUE(reader.has_topology());
    EXPECT_TRUE(reader.has_trajectory());
    EXPECT_FALSE(reader.has_trajectory_metadata());
    EXPECT_TRUE(reader.has_bonds());
    EXPECT_TRUE(reader.can_read(".mol2"));
    EXPECT_FALSE(reader.can_read(".psf"));
    EXPECT_FALSE(reader.can_read(".mol"));
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
    EXPECT_EQ(Atom(0, 0, atoms).name(), "N");
    EXPECT_EQ(Atom(1, 0, atoms).name(), "CA");

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
}

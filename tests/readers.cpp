#include "matchers.hpp"
#include "readers/MolReader.hpp"
#include "readers/MolfileReader.hpp"
#include "core/MolData.hpp"
#include <molpp/MolError.hpp>
#include <molpp/Atom.hpp>
#include <molpp/AtomSel.hpp>
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <optional>

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
    EXPECT_FALSE(MolfileReader::can_read(""));
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
    for (index_t i = 0; i < data->size(); ++i)
    {
        atoms.push_back(Atom(i, 0, data));
    }
    std::vector<mol::Atom> m2_atoms;
    for (index_t i = 0; i < m2_data->size(); ++i)
    {
        m2_atoms.push_back(Atom(i, std::nullopt, m2_data));
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
    ASSERT_EQ(data->trajectory().num_frames(), 0);

    ASSERT_EQ(reader.read_timestep(data), MolReader::SUCCESS);
    ASSERT_EQ(data->trajectory().num_frames(), 1);
    EXPECT_THAT(data->trajectory().timestep(0).coords().reshaped(), ElementsAre(2, -2, 0, 4, -4, 2));

    ASSERT_EQ(reader.skip_timestep(data), MolReader::SUCCESS);
    ASSERT_EQ(data->trajectory().num_frames(), 1);

    ASSERT_EQ(reader.read_timestep(data), MolReader::SUCCESS);
    ASSERT_EQ(data->trajectory().num_frames(), 2);
    EXPECT_THAT(data->trajectory().timestep(1).coords().reshaped(), ElementsAre(24, -24, 0, 48, -48, 24));

    ASSERT_EQ(reader.read_timestep(data), MolReader::END);
    ASSERT_EQ(reader.skip_timestep(data), MolReader::END);

    reader.close();

    /*
     * Bonds and charges
     */
    // Bond partners
    EXPECT_THAT(AtomSel(m2_atoms[0]).bonded().indices(), ElementsAre(0, 1));
    EXPECT_THAT(AtomSel(m2_atoms[1]).bonded().indices(), ElementsAre(0, 1, 2, 3));
    EXPECT_THAT(AtomSel(m2_atoms[2]).bonded().indices(), ElementsAre(1, 2, 4, 7));
    EXPECT_THAT(AtomSel(m2_atoms[3]).bonded().indices(), ElementsAre(1, 3, 5, 8));
    EXPECT_THAT(AtomSel(m2_atoms[4]).bonded().indices(), ElementsAre(2, 4, 6, 9));
    EXPECT_THAT(AtomSel(m2_atoms[5]).bonded().indices(), ElementsAre(3, 5, 6, 10));
    EXPECT_THAT(AtomSel(m2_atoms[6]).bonded().indices(), ElementsAre(4, 5, 6, 11));
    EXPECT_THAT(AtomSel(m2_atoms[7]).bonded().indices(), ElementsAre(2, 7));
    EXPECT_THAT(AtomSel(m2_atoms[8]).bonded().indices(), ElementsAre(3, 8));
    EXPECT_THAT(AtomSel(m2_atoms[9]).bonded().indices(), ElementsAre(4, 9));
    EXPECT_THAT(AtomSel(m2_atoms[10]).bonded().indices(), ElementsAre(5, 10));
    EXPECT_THAT(AtomSel(m2_atoms[11]).bonded().indices(), ElementsAre(6, 11));

    // Bond orders
    for (Atom &atom : m2_atoms)
    {
        for (auto bond : atom.bonds())
        {
            EXPECT_FALSE(bond->guessed());
        }
    }
    EXPECT_EQ(m2_atoms[0].bond(1)->order(), 1);
    EXPECT_FALSE(m2_atoms[0].bond(1)->guessed_order());
    EXPECT_EQ(m2_atoms[1].bond(2)->order(), 0);
    EXPECT_TRUE(m2_atoms[1].bond(2)->guessed_order());
    EXPECT_EQ(m2_atoms[1].bond(3)->order(), 0);
    EXPECT_TRUE(m2_atoms[1].bond(3)->guessed_order());
    EXPECT_EQ(m2_atoms[2].bond(4)->order(), 0);
    EXPECT_TRUE(m2_atoms[2].bond(4)->guessed_order());
    EXPECT_EQ(m2_atoms[2].bond(7)->order(), 1);
    EXPECT_FALSE(m2_atoms[2].bond(7)->guessed_order());
    EXPECT_EQ(m2_atoms[3].bond(5)->order(), 0);
    EXPECT_TRUE(m2_atoms[3].bond(5)->guessed_order());
    EXPECT_EQ(m2_atoms[3].bond(8)->order(), 1);
    EXPECT_FALSE(m2_atoms[3].bond(8)->guessed_order());
    EXPECT_EQ(m2_atoms[4].bond(6)->order(), 0);
    EXPECT_TRUE(m2_atoms[4].bond(6)->guessed_order());
    EXPECT_EQ(m2_atoms[4].bond(9)->order(), 1);
    EXPECT_FALSE(m2_atoms[4].bond(9)->guessed_order());
    EXPECT_EQ(m2_atoms[5].bond(6)->order(), 0);
    EXPECT_TRUE(m2_atoms[5].bond(6)->guessed_order());
    EXPECT_EQ(m2_atoms[5].bond(10)->order(), 1);
    EXPECT_FALSE(m2_atoms[5].bond(10)->guessed_order());
    EXPECT_EQ(m2_atoms[6].bond(11)->order(), 1);
    EXPECT_FALSE(m2_atoms[6].bond(11)->guessed_order());

    /*
     * Residue detection
     */
    reader = MolfileReader(".pdb");
    ASSERT_EQ(reader.open("4lad.pdb"), MolReader::SUCCESS);
    data = reader.read_atoms();
    reader.close();
    ASSERT_THAT(data, NotNull());
    EXPECT_EQ(data->atoms().residue(0), 0);
    EXPECT_EQ(data->atoms().residue(3), 0);
    EXPECT_EQ(data->atoms().residue(4), 1);
    EXPECT_EQ(data->atoms().residue(1230), 153);
    EXPECT_EQ(data->atoms().residue(1231), 154);
    EXPECT_EQ(data->atoms().residue(1561), 195);
    EXPECT_EQ(data->atoms().residue(1562), 196);
    EXPECT_EQ(data->atoms().residue(1790), 222);
    EXPECT_EQ(data->atoms().residue(1791), 223);
    EXPECT_EQ(data->atoms().residue(1792), 224);
    EXPECT_EQ(data->atoms().residue(1793), 225);
    EXPECT_EQ(data->atoms().residue(1798), 225);
    EXPECT_EQ(data->atoms().residue(1799), 226);
    EXPECT_EQ(data->atoms().residue(1804), 226);
    // Water molecules
    for (index_t i = 0; i < 39; ++i)
    {
        EXPECT_EQ(data->atoms().residue(1805 + i), 227 + i);
    }
}

TEST(Readers, PDB) {
    ASSERT_TRUE(MolfileReader::can_read(".pdb"));
    ASSERT_TRUE(MolfileReader::can_read(".ent"));
    MolfileReader reader(".pdb");
    EXPECT_TRUE(reader.has_topology());
    EXPECT_TRUE(reader.has_trajectory());
    EXPECT_FALSE(reader.has_trajectory_metadata());
    EXPECT_TRUE(reader.has_bonds());
}

TEST(Readers, Mol2) {
    ASSERT_TRUE(MolfileReader::can_read(".mol2"));
    MolfileReader reader(".mol2");
    EXPECT_TRUE(reader.has_topology());
    EXPECT_TRUE(reader.has_trajectory());
    EXPECT_FALSE(reader.has_trajectory_metadata());
    EXPECT_TRUE(reader.has_bonds());
}

TEST(Readers, PSF) {
    ASSERT_TRUE(MolfileReader::can_read(".psf"));
    MolfileReader reader(".psf");
    EXPECT_TRUE(reader.has_topology());
    EXPECT_FALSE(reader.has_trajectory());
    EXPECT_FALSE(reader.has_trajectory_metadata());
    EXPECT_TRUE(reader.has_bonds());

    // Read topology
    ASSERT_EQ(reader.open("dipeptide.psf"), MolReader::SUCCESS);
    auto data = reader.read_atoms();
    reader.close();
    ASSERT_THAT(data, NotNull());
    ASSERT_EQ(data->size(), 22);

    // Atom properties
    std::vector<mol::Atom> atoms;
    for (index_t i = 0; i < data->size(); ++i)
    {
        atoms.push_back(Atom(i, std::nullopt, data));
    }
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::resid), {129, 129, 129, 129, 130, 130, 130, 130, 130, 130, 130, 129, 129, 129, 129, 130, 130, 130, 130, 130, 130, 130}));
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::residue_id), {0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3}));
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::atomic), {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}));
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::occupancy), {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}));
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::tempfactor), {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}));
    EXPECT_THAT(atoms, Pointwise(PropFloat(&Atom::mass, 1e-5), {14.007, 12.011, 12.011, 15.999, 14.007, 12.011, 12.011, 12.011, 12.011, 12.011, 15.999, 14.007, 12.011, 12.011, 15.999, 14.007, 12.011, 12.011, 12.011, 12.011, 12.011, 15.999}));
    EXPECT_THAT(atoms, Pointwise(PropFloat(&Atom::radius, 1e-5), {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}));
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::name), {"N", "CA", "C", "O", "N", "CD", "CA", "CB", "CG", "C", "O", "N", "CA", "C", "O", "N", "CD", "CA", "CB", "CG", "C", "O"}));
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::type), {"NH1", "CT2", "C", "O", "N", "CP3", "CP1", "CP2", "CP2", "C", "O", "NH1", "CT2", "C", "O", "N", "CP3", "CP1", "CP2", "CP2", "C", "O"}));
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::resname), {"GLY", "GLY", "GLY", "GLY", "PRO", "PRO", "PRO", "PRO", "PRO", "PRO", "PRO", "GLY", "GLY", "GLY", "GLY", "PRO", "PRO", "PRO", "PRO", "PRO", "PRO", "PRO"}));
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::segid), {"KCH1", "KCH1", "KCH1", "KCH1", "KCH1", "KCH1", "KCH1", "KCH1", "KCH1", "KCH1", "KCH1", "KCH3", "KCH3", "KCH3", "KCH3", "KCH3", "KCH3", "KCH3", "KCH3", "KCH3", "KCH3", "KCH3"}));
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::chain), {"K", "K", "K", "K", "K", "K", "K", "K", "K", "K", "K", "K", "K", "K", "K", "K", "K", "K", "K", "K", "K", "K"}));
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::altloc), {"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""}));
    EXPECT_THAT(atoms, Pointwise(PropFloat(&Atom::charge, 1e-5), {-0.47, -0.02, 0.51, -0.51, -0.29, 0.0, 0.02, -0.18, -0.18, 0.51, -0.51, -0.47, -0.02, 0.51, -0.51, -0.29, 0.0, 0.02, -0.18, -0.18, 0.51, -0.51}));

    // Bonds
    EXPECT_THAT(AtomSel(atoms[0]).bonded().indices(), ElementsAre(0, 1));
    EXPECT_THAT(AtomSel(atoms[1]).bonded().indices(), ElementsAre(0, 1, 2));
    EXPECT_THAT(AtomSel(atoms[2]).bonded().indices(), ElementsAre(1, 2, 3, 4));
    EXPECT_THAT(AtomSel(atoms[3]).bonded().indices(), ElementsAre(2, 3));
    EXPECT_THAT(AtomSel(atoms[4]).bonded().indices(), ElementsAre(2, 4, 5, 6));
    EXPECT_THAT(AtomSel(atoms[5]).bonded().indices(), ElementsAre(4, 5, 8));
    EXPECT_THAT(AtomSel(atoms[6]).bonded().indices(), ElementsAre(4, 6, 7, 9));
    EXPECT_THAT(AtomSel(atoms[7]).bonded().indices(), ElementsAre(6, 7, 8));
    EXPECT_THAT(AtomSel(atoms[8]).bonded().indices(), ElementsAre(5, 7, 8));
    EXPECT_THAT(AtomSel(atoms[9]).bonded().indices(), ElementsAre(6, 9, 10));
    EXPECT_THAT(AtomSel(atoms[10]).bonded().indices(), ElementsAre(9, 10));
    EXPECT_THAT(AtomSel(atoms[11]).bonded().indices(), ElementsAre(11, 12));
    EXPECT_THAT(AtomSel(atoms[12]).bonded().indices(), ElementsAre(11, 12, 13));
    EXPECT_THAT(AtomSel(atoms[13]).bonded().indices(), ElementsAre(12, 13, 14, 15));
    EXPECT_THAT(AtomSel(atoms[14]).bonded().indices(), ElementsAre(13, 14));
    EXPECT_THAT(AtomSel(atoms[15]).bonded().indices(), ElementsAre(13, 15, 16, 17));
    EXPECT_THAT(AtomSel(atoms[16]).bonded().indices(), ElementsAre(15, 16, 19));
    EXPECT_THAT(AtomSel(atoms[17]).bonded().indices(), ElementsAre(15, 17, 18, 20));
    EXPECT_THAT(AtomSel(atoms[18]).bonded().indices(), ElementsAre(17, 18, 19));
    EXPECT_THAT(AtomSel(atoms[19]).bonded().indices(), ElementsAre(16, 18, 19));
    EXPECT_THAT(AtomSel(atoms[20]).bonded().indices(), ElementsAre(17, 20, 21));
    EXPECT_THAT(AtomSel(atoms[21]).bonded().indices(), ElementsAre(20, 21));
}

TEST(Readers, XTC) {
    ASSERT_TRUE(MolfileReader::can_read(".psf"));
    ASSERT_TRUE(MolfileReader::can_read(".xtc"));
    MolfileReader xtc_reader(".xtc");
    EXPECT_FALSE(xtc_reader.has_topology());
    EXPECT_TRUE(xtc_reader.has_trajectory());
    EXPECT_FALSE(xtc_reader.has_trajectory_metadata());
    EXPECT_FALSE(xtc_reader.has_bonds());

    // Load the topology
    MolfileReader psf_reader(".psf");
    ASSERT_EQ(psf_reader.open("dipeptide.psf"), MolReader::SUCCESS);
    auto data = psf_reader.read_atoms();
    psf_reader.close();
    ASSERT_THAT(data, NotNull());
    ASSERT_EQ(data->size(), 22);

    // Load the trajectory
    ASSERT_EQ(xtc_reader.open("dipeptide.xtc"), MolReader::SUCCESS);
    ASSERT_EQ(xtc_reader.check_timestep_read(data), MolReader::SUCCESS);
    ASSERT_EQ(xtc_reader.read_timestep(data), MolReader::SUCCESS);
    ASSERT_EQ(xtc_reader.read_timestep(data), MolReader::SUCCESS);
    ASSERT_EQ(data->trajectory().num_frames(), 2);
    xtc_reader.close();

    // Check coordinates
    EXPECT_THAT(data->trajectory().timestep(0).coords().reshaped(), Pointwise(FloatNear(1e-5), {488.93002, 780.79004, -2.42000, 488.31000, 779.74005, -1.57000, 488.95001, 778.40002, -1.68000, 489.36002, 778.07007, -2.80000, 488.83002, 777.53003, -0.60000, 488.71002, 777.94006, 0.79000, 489.25003, 776.12006, -0.75000, 489.37006, 775.66003, 0.70000, 489.78000, 776.99005, 1.34000, 488.46002, 775.29004, -1.63000, 487.31003, 775.13000, -1.27000, 438.54004, 830.56006, -3.76000, 439.46002, 831.70007, -3.88000, 440.85004, 831.40009, -4.47000, 441.05005, 830.21002, -4.75000, 441.64001, 832.37006, -4.79000, 441.41003, 833.80005, -4.52000, 443.00003, 832.08002, -5.41000, 443.75003, 833.43005, -5.27000, 442.56000, 834.47009, -5.24000, 443.85004, 831.02002, -4.66000, 443.70001, 830.90002, -3.40000}));
    EXPECT_THAT(data->trajectory().timestep(1).coords().reshaped(), Pointwise(FloatNear(1e-5), {488.76004, 782.39008, -1.81000, 487.96002, 781.47003, -1.06000, 488.48004, 780.12000, -1.26000, 488.72003, 779.78003, -2.39000, 488.78003, 779.27002, -0.28000, 488.90002, 779.65002, 1.13000, 489.07001, 777.84003, -0.42000, 489.22000, 777.34998, 0.98000, 489.68002, 778.56006, 1.77000, 488.13004, 776.95007, -1.27000, 486.92001, 776.92004, -1.05000, 437.76001, 830.54004, -4.03000, 438.49002, 831.81006, -3.87000, 439.61002, 832.06000, -4.89000, 439.85001, 831.20001, -5.72000, 440.36005, 833.23004, -4.93000, 439.95001, 834.44000, -4.13000, 441.73001, 833.40002, -5.36000, 442.09003, 834.89008, -5.01000, 440.76001, 835.60004, -4.67000, 442.71002, 832.34998, -4.79000, 442.83002, 832.17004, -3.55000}));
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
    EXPECT_EQ(pdb_reader->read_trajectory("", atoms), MolReader::FAILED);
    EXPECT_EQ(pdb_reader->read_trajectory("tiny.pdb", atoms), MolReader::WRONG_ATOMS);

    EXPECT_EQ(pdb_reader->read_trajectory("traj.pdb", atoms), MolReader::SUCCESS);
    EXPECT_EQ(atoms->trajectory().num_frames(), 4);
    EXPECT_EQ(pdb_reader->read_trajectory("traj.pdb", atoms, 0, 1, 0), MolReader::SUCCESS);
    EXPECT_EQ(atoms->trajectory().num_frames(), 5);
    EXPECT_EQ(pdb_reader->read_trajectory("traj.pdb", atoms, 1, 2, 2), MolReader::SUCCESS);
    EXPECT_EQ(atoms->trajectory().num_frames(), 6);
    EXPECT_EQ(pdb_reader->read_trajectory("traj.pdb", atoms, 2, 0), MolReader::SUCCESS);
    EXPECT_EQ(pdb_reader->read_trajectory("traj.pdb", atoms, -2, 0, 2), MolReader::SUCCESS);
    EXPECT_EQ(atoms->trajectory().num_frames(), 6);
}

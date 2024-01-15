#include "readers/MolfileReader.hpp"
#include <molpp/Atom.hpp>
#include <molpp/Residue.hpp>
#include <molpp/internal/MolData.hpp>
#include <molpp/AtomSel.hpp>
#include <molpp/Property.hpp>
#include <molpp/MolError.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace testing;
using namespace mol;
using namespace mol::internal;

class BaseMolfileReaderTest : public ::testing::Test
{
protected:
    static std::unique_ptr<MolData> read(std::string const& filename, std::string const& extension)
    {
        MolfileReader reader(extension);
        reader = MolfileReader(extension);
        reader.open(filename);
        auto data = reader.read_atoms();
        reader.close();
        return std::move(data);
    }
};

TEST(MolfileReaderTest, UnknownExtension)
{
    EXPECT_THROW(MolfileReader(""), MolError);
    EXPECT_THROW(MolfileReader("bar"), MolError);
    EXPECT_THROW(MolfileReader(".foo"), MolError);
}

TEST(MolfileReaderTest, InvalidCanRead)
{
    EXPECT_FALSE(MolfileReader::can_read(""));
}

TEST(MolfileReaderTest, ReadBeforeOpen)
{
    MolfileReader reader(".pdb");
    EXPECT_THROW(reader.read_atoms(), MolError);
}

TEST(MolfileReaderTest, OpenValid)
{
    MolfileReader reader(".pdb");
    EXPECT_EQ(reader.open("tiny.pdb"), MolReader::SUCCESS);
    reader.close();
}

TEST(MolfileReaderTest, OpenTwice)
{
    MolfileReader reader(".pdb");
    ASSERT_EQ(reader.open("tiny.pdb"), MolReader::SUCCESS);
    EXPECT_EQ(reader.open("tiny.pdb"), MolReader::INVALID);
    reader.close();
}

TEST(MolfileReaderTest, OpenInvalid)
{
    MolfileReader reader(".pdb");
    EXPECT_EQ(reader.open(""), MolReader::FAILED);
}

struct FormatInfo
{
    std::string const format;
    bool const has_topology;
    bool const has_trajectory;
    bool const has_trajectory_metadata;
    bool const has_bonds;
};

class MolfileReaderFormatTest : public testing::TestWithParam<FormatInfo>
{
};

TEST_P(MolfileReaderFormatTest, CanRead)
{
    EXPECT_TRUE(MolfileReader::can_read(GetParam().format));
}

TEST_P(MolfileReaderFormatTest, ConstructReader)
{
    EXPECT_NO_THROW(MolfileReader(GetParam().format));
}

TEST_P(MolfileReaderFormatTest, HasTopology)
{
    MolfileReader reader(GetParam().format);
    EXPECT_EQ(reader.has_topology(), GetParam().has_topology);
}

TEST_P(MolfileReaderFormatTest, HasTrajectory)
{
    MolfileReader reader(GetParam().format);
    EXPECT_EQ(reader.has_trajectory(), GetParam().has_trajectory);
}

TEST_P(MolfileReaderFormatTest, HasTrajectoryMetadata)
{
    MolfileReader reader(GetParam().format);
    EXPECT_EQ(reader.has_trajectory_metadata(), GetParam().has_trajectory_metadata);
}

TEST_P(MolfileReaderFormatTest, HasBonds)
{
    MolfileReader reader(GetParam().format);
    EXPECT_EQ(reader.has_bonds(), GetParam().has_bonds);
}

INSTANTIATE_TEST_SUITE_P(Formats, MolfileReaderFormatTest, Values(
    FormatInfo{".pdb", true, true, false, true},
    FormatInfo{".ent", true, true, false, true},
    FormatInfo{".mol2", true, true, false, true},
    FormatInfo{".psf", true, false, false, true},
    FormatInfo{".xtc", false, true, false, false}
));

struct TopologyInfo
{
    std::string const filename;
    std::string const extension;
    size_t const num_atoms;
};

class MolfileReaderTopologyReadingTest : public testing::TestWithParam<TopologyInfo>
{
};

TEST_P(MolfileReaderTopologyReadingTest, Read)
{
    MolfileReader reader(GetParam().extension);
    ASSERT_EQ(reader.open(GetParam().filename), MolReader::SUCCESS);
    auto data = reader.read_atoms();
    reader.close();
    ASSERT_THAT(data, NotNull());
    EXPECT_EQ(data->properties().size<Atom>(), GetParam().num_atoms);
}

INSTANTIATE_TEST_SUITE_P(Files, MolfileReaderTopologyReadingTest, Values(
    TopologyInfo{"tiny.pdb", ".pdb", 6},
    TopologyInfo{"fluorobenzene.mol2", ".mol2", 12},
    TopologyInfo{"dipeptide.psf", ".psf", 22}
));

class MolfileReaderTrajectoryFunctionsTest : public BaseMolfileReaderTest
{
protected:
    void SetUp() override
    {
        data = read("tiny.pdb", ".pdb");
        ASSERT_THAT(data, NotNull());
        ASSERT_EQ(data->properties().size<Atom>(), 6);
    }

    std::unique_ptr<MolData> data;
};

TEST_F(MolfileReaderTrajectoryFunctionsTest, WrongAtomCount)
{
    MolfileReader reader(".pdb");
    ASSERT_EQ(reader.open("traj.pdb"), MolReader::SUCCESS);
    EXPECT_EQ(reader.check_timestep_read(*data), MolReader::WRONG_ATOMS);
}

TEST_F(MolfileReaderTrajectoryFunctionsTest, ClosedReader)
{
    MolfileReader reader(".pdb");
    EXPECT_EQ(reader.check_timestep_read(*data), MolReader::INVALID);
}

TEST_F(MolfileReaderTrajectoryFunctionsTest, SkipTimestep)
{
    MolfileReader reader(".pdb");
    reader.open("traj.pdb");
    EXPECT_EQ(reader.skip_timestep(*data), MolReader::SUCCESS);
    EXPECT_EQ(data->properties().num_frames(), 0);
    reader.close();
}

TEST_F(MolfileReaderTrajectoryFunctionsTest, ReadTimestep)
{
    MolfileReader reader(".pdb");
    reader.open("traj.pdb");
    EXPECT_EQ(reader.read_timestep(*data), MolReader::SUCCESS);
    EXPECT_EQ(data->properties().num_frames(), 1);
    reader.close();
}

TEST_F(MolfileReaderTrajectoryFunctionsTest, ReadAndSkipTimestep)
{
    MolfileReader reader(".pdb");
    reader.open("traj.pdb");
    reader.skip_timestep(*data);
    reader.read_timestep(*data);
    reader.skip_timestep(*data);
    reader.read_timestep(*data);
    EXPECT_EQ(data->properties().num_frames(), 2);
    reader.close();
}

TEST_F(MolfileReaderTrajectoryFunctionsTest, ReadPastEnd)
{
    MolfileReader reader(".pdb");
    reader.open("traj.pdb");
    reader.skip_timestep(*data);
    reader.skip_timestep(*data);
    reader.skip_timestep(*data);
    reader.skip_timestep(*data);
    EXPECT_EQ(reader.read_timestep(*data), MolReader::END);
    reader.close();
}

TEST_F(MolfileReaderTrajectoryFunctionsTest, SkipPastEnd)
{
    MolfileReader reader(".pdb");
    reader.open("traj.pdb");
    reader.skip_timestep(*data);
    reader.skip_timestep(*data);
    reader.skip_timestep(*data);
    reader.skip_timestep(*data);
    EXPECT_EQ(reader.skip_timestep(*data), MolReader::END);
    reader.close();
}

struct TrajectoryInfo
{
    std::string const topology;
    std::string const topology_ext;
    std::string const trajectory;
    std::string const trajectory_ext;
    size_t const num_frames;
    std::vector<std::vector<double>> const positions;
};

class MolfileReaderReadTrajectoryTest : public testing::TestWithParam<TrajectoryInfo>
{
};

TEST_P(MolfileReaderReadTrajectoryTest, Read)
{
    MolfileReader reader(GetParam().topology_ext);
    ASSERT_EQ(reader.open(GetParam().topology), MolReader::SUCCESS);
    auto data = reader.read_atoms();
    reader.close();
    ASSERT_THAT(data, NotNull());

    ASSERT_NO_THROW(reader = MolfileReader(GetParam().trajectory_ext));
    ASSERT_EQ(reader.open(GetParam().trajectory), MolReader::SUCCESS);
    ASSERT_EQ(reader.check_timestep_read(*data), MolReader::SUCCESS);
    for (size_t i = 0; i < GetParam().num_frames; i++)
    {
        EXPECT_EQ(reader.read_timestep(*data), MolReader::SUCCESS);
    }
    reader.close();

    EXPECT_EQ(data->properties().num_frames(), GetParam().num_frames);

    auto expected_positions = GetParam().positions;
    for (size_t frame = 0; frame < GetParam().num_frames; frame++)
    {
        Position* positions = data->properties().get<Atom, Position>(frame);
        ASSERT_THAT(positions, NotNull());
        EXPECT_THAT(positions->positions().reshaped(), Pointwise(FloatNear(1e-5), expected_positions[frame])) << "Frame " << frame;
    }
}

INSTANTIATE_TEST_SUITE_P(Files, MolfileReaderReadTrajectoryTest, Values(
    TrajectoryInfo{"traj.pdb", ".pdb", "traj.pdb", ".pdb", 4, {{1, -1, 0, 2, -2, 1}, {2, -2, 0, 4, -4, 2}, {6, -6, 0, 12, -12, 6}, {24, -24, 0, 48, -48, 24}}},
    TrajectoryInfo{"fluorobenzene.mol2", ".mol2", "fluorobenzene.mol2", ".mol2", 1, {{-2.3436,  0.0000, -.0001, -1.0043,  0.0000, 0.0002, -0.30690,  1.2080,  0.0000, -0.30680, -1.2079,  0.0000,  1.0880,  1.2079,  0.0000, 1.0880, -1.2080,  0.0000,  1.7855,  0.0000, 0.0000, -0.85010,  2.1483, -0.0001, -0.85010, -2.1483, -0.0001,  1.6310,  2.1485, -0.0001, 1.6311, -2.1485,  0.0000,  2.8716,  0.0001, 0.0000}}},
    TrajectoryInfo{"dipeptide.psf", ".psf", "dipeptide.xtc", ".xtc", 2, {{488.93002, 780.79004, -2.42000, 488.31000, 779.74005, -1.57000, 488.95001, 778.40002, -1.68000, 489.36002, 778.07007, -2.80000, 488.83002, 777.53003, -0.60000, 488.71002, 777.94006, 0.79000, 489.25003, 776.12006, -0.75000, 489.37006, 775.66003, 0.70000, 489.78000, 776.99005, 1.34000, 488.46002, 775.29004, -1.63000, 487.31003, 775.13000, -1.27000, 438.54004, 830.56006, -3.76000, 439.46002, 831.70007, -3.88000, 440.85004, 831.40009, -4.47000, 441.05005, 830.21002, -4.75000, 441.64001, 832.37006, -4.79000, 441.41003, 833.80005, -4.52000, 443.00003, 832.08002, -5.41000, 443.75003, 833.43005, -5.27000, 442.56000, 834.47009, -5.24000, 443.85004, 831.02002, -4.66000, 443.70001, 830.90002, -3.40000}, {488.76004, 782.39008, -1.81000, 487.96002, 781.47003, -1.06000, 488.48004, 780.12000, -1.26000, 488.72003, 779.78003, -2.39000, 488.78003, 779.27002, -0.28000, 488.90002, 779.65002, 1.13000, 489.07001, 777.84003, -0.42000, 489.22000, 777.34998, 0.98000, 489.68002, 778.56006, 1.77000, 488.13004, 776.95007, -1.27000, 486.92001, 776.92004, -1.05000, 437.76001, 830.54004, -4.03000, 438.49002, 831.81006, -3.87000, 439.61002, 832.06000, -4.89000, 439.85001, 831.20001, -5.72000, 440.36005, 833.23004, -4.93000, 439.95001, 834.44000, -4.13000, 441.73001, 833.40002, -5.36000, 442.09003, 834.89008, -5.01000, 440.76001, 835.60004, -4.67000, 442.71002, 832.34998, -4.79000, 442.83002, 832.17004, -3.55000}}}
));

template<IsAtomAggregate Aggregate, IsProperty PropertyType, size_t Size>
void test_property(std::array<typename PropertyType::value_type, Size> const expected, MolData* data)
{
    for (size_t i = 0; i < Size; i++)
    {
        Aggregate aggregate(i, std::nullopt, data);
        EXPECT_EQ(aggregate.template get<PropertyType>(), expected[i]) << i;
    }
}

template<IsAtomAggregate Aggregate, IsProperty PropertyType>
void test_no_property(MolData* data)
{
    size_t num_elements = data->properties().size<Aggregate>();
    for (size_t i = 0; i < num_elements; i++)
    {
        Aggregate aggregate(i, std::nullopt, data);
        EXPECT_FALSE(aggregate.template has<PropertyType>()) << i;
    }
}

class PDBMolfileReaderTest : public BaseMolfileReaderTest
{
protected:
    void SetUp() override
    {
        ASSERT_THAT(pdb(), NotNull());
        ASSERT_EQ(pdb()->properties().size<Atom>(), 6);
    }

    MolData* pdb()
    {
        static std::unique_ptr<MolData> data = read("tiny.pdb", ".pdb");
        return data.get();
    }
};

TEST_F(PDBMolfileReaderTest, AtomName)
{
    test_property<Atom, Name, 6>({"N", "NA", "O1", "O22", "ND", "ND"}, pdb());
}

TEST_F(PDBMolfileReaderTest, AtomType)
{
    test_property<Atom, Type, 6>({"N", "NA", "O1", "O22", "ND", "ND"}, pdb());
}

TEST_F(PDBMolfileReaderTest, Occupancy)
{
    test_property<Atom, Occupancy, 6>({999, 1, 1, 1, 1, 1}, pdb());
}

TEST_F(PDBMolfileReaderTest, TemperatureFactor)
{
    test_property<Atom, TemperatureFactor, 6>({-99, -1, -1, -1, -1, -1}, pdb());
}

TEST_F(PDBMolfileReaderTest, Mass)
{
    test_property<Atom, Mass, 6>({14.0067, 14.0067, 15.9994, 15.9994, 14.0067, 14.0067}, pdb());
}

TEST_F(PDBMolfileReaderTest, Radius)
{
    test_property<Atom, Radius, 6>({1.55, 1.55, 1.52, 1.52, 1.55, 1.55}, pdb());
}

TEST_F(PDBMolfileReaderTest, AtomicNumber)
{
    test_property<Atom, AtomicNumber, 6>({7, 7, 8, 8, 7, 7}, pdb());
}

TEST_F(PDBMolfileReaderTest, Charge)
{
    test_no_property<Atom, Charge>(pdb());
}

TEST_F(PDBMolfileReaderTest, AlternateLocation)
{
    test_property<Atom, AlternateLocation, 6>({" ", " ", " ", " ", "A", "B"}, pdb());
}

TEST_F(PDBMolfileReaderTest, InsertionCode)
{
    test_property<Atom, InsertionCode, 6>({"A", "B", " ", " ", " ", " "}, pdb());
}

TEST_F(PDBMolfileReaderTest, ResID)
{
    test_property<Residue, ResID, 5>({3, 339, 201, 801, 85}, pdb());
}

TEST_F(PDBMolfileReaderTest, ResName)
{
    test_property<Residue, ResName, 5>({"GLY", "ASP", "HOH", "HOH", "ASP"}, pdb());
}

TEST_F(PDBMolfileReaderTest, Bonded)
{
    EXPECT_THAT(AtomSel(Atom(0, std::nullopt, pdb())).bonded().indices(), ElementsAre(0, 2, 3));
    EXPECT_THAT(AtomSel(Atom(2, std::nullopt, pdb())).bonded().indices(), ElementsAre(0, 2));
    EXPECT_THAT(AtomSel(Atom(3, std::nullopt, pdb())).bonded().indices(), ElementsAre(0, 3));
}

TEST_F(PDBMolfileReaderTest, GuessedBond)
{
    size_t const num_atoms = pdb()->properties().size<Atom>();
    for (size_t i = 0; i < num_atoms; i++)
    {
        for (auto bond : Atom(i, std::nullopt, pdb()).bonds())
        {
            EXPECT_FALSE(bond->guessed());
        }
    }
}

TEST_F(PDBMolfileReaderTest, BondOrder)
{
    size_t const num_atoms = pdb()->properties().size<Atom>();
    for (size_t i = 0; i < num_atoms; i++)
    {
        for (auto bond : Atom(i, std::nullopt, pdb()).bonds())
        {
            EXPECT_EQ(bond->order(), 1) << bond->atom1() << "-" << bond->atom2();
        }
    }
}

TEST_F(PDBMolfileReaderTest, GuessedOrder)
{
    size_t const num_atoms = pdb()->properties().size<Atom>();
    for (size_t i = 0; i < num_atoms; i++)
    {
        for (auto bond : Atom(i, std::nullopt, pdb()).bonds())
        {
            EXPECT_TRUE(bond->guessed_order()) << bond->atom1() << "-" << bond->atom2();
        }
    }
}

class Mol2MolfileReaderTest : public BaseMolfileReaderTest
{
protected:
    void SetUp() override
    {
        ASSERT_THAT(mol2(), NotNull());
        ASSERT_EQ(mol2()->properties().size<Atom>(), 12);
    }

    MolData* mol2()
    {
        static std::unique_ptr<MolData> data = read("fluorobenzene.mol2", ".mol2");
        return data.get();
    }
};

TEST_F(Mol2MolfileReaderTest, AtomName)
{
    test_property<Atom, Name, 12>({"F", "C", "C", "C", "C", "C", "C", "H", "H", "H", "H", "H"}, mol2());
}

TEST_F(Mol2MolfileReaderTest, AtomType)
{
    test_property<Atom, Type, 12>({"F", "C.ar", "C.ar", "C.ar", "C.ar", "C.ar", "C.ar", "H", "H", "H", "H", "H"}, mol2());
}

TEST_F(Mol2MolfileReaderTest, Occupancy)
{
    test_no_property<Atom, Occupancy>(mol2());
}

TEST_F(Mol2MolfileReaderTest, TemperatureFactor)
{
    test_no_property<Atom, TemperatureFactor>(mol2());
}

TEST_F(Mol2MolfileReaderTest, Mass)
{
    test_no_property<Atom, Mass>(mol2());
}

TEST_F(Mol2MolfileReaderTest, Radius)
{
    test_no_property<Atom, Radius>(mol2());
}

TEST_F(Mol2MolfileReaderTest, AtomicNumber)
{
    test_no_property<Atom, AtomicNumber>(mol2());
}

TEST_F(Mol2MolfileReaderTest, Charge)
{
    test_property<Atom, Charge, 12>({-0.2055,  0.1234, -0.0265, -0.0265, -0.0590, -0.0590, -0.0616,  0.0646, 0.0646,  0.0618,  0.0618,  0.0618}, mol2());
}

TEST_F(Mol2MolfileReaderTest, AlternateLocation)
{
    test_no_property<Atom, AlternateLocation>(mol2());
}

TEST_F(Mol2MolfileReaderTest, InsertionCode)
{
    test_no_property<Atom, InsertionCode>(mol2());
}

TEST_F(Mol2MolfileReaderTest, ResID)
{
    test_property<Residue, ResID, 1>({1}, mol2());
}

TEST_F(Mol2MolfileReaderTest, ResName)
{
    test_property<Residue, ResName, 1>({"FLB1"}, mol2());
}

TEST_F(Mol2MolfileReaderTest, Bonded)
{
    EXPECT_THAT(AtomSel(Atom(0, std::nullopt, mol2())).bonded().indices(), ElementsAre(0, 1));
    EXPECT_THAT(AtomSel(Atom(1, std::nullopt, mol2())).bonded().indices(), ElementsAre(0, 1, 2, 3));
    EXPECT_THAT(AtomSel(Atom(2, std::nullopt, mol2())).bonded().indices(), ElementsAre(1, 2, 4, 7));
    EXPECT_THAT(AtomSel(Atom(3, std::nullopt, mol2())).bonded().indices(), ElementsAre(1, 3, 5, 8));
    EXPECT_THAT(AtomSel(Atom(4, std::nullopt, mol2())).bonded().indices(), ElementsAre(2, 4, 6, 9));
    EXPECT_THAT(AtomSel(Atom(5, std::nullopt, mol2())).bonded().indices(), ElementsAre(3, 5, 6, 10));
    EXPECT_THAT(AtomSel(Atom(6, std::nullopt, mol2())).bonded().indices(), ElementsAre(4, 5, 6, 11));
    EXPECT_THAT(AtomSel(Atom(7, std::nullopt, mol2())).bonded().indices(), ElementsAre(2, 7));
    EXPECT_THAT(AtomSel(Atom(8, std::nullopt, mol2())).bonded().indices(), ElementsAre(3, 8));
    EXPECT_THAT(AtomSel(Atom(9, std::nullopt, mol2())).bonded().indices(), ElementsAre(4, 9));
    EXPECT_THAT(AtomSel(Atom(10, std::nullopt, mol2())).bonded().indices(), ElementsAre(5, 10));
    EXPECT_THAT(AtomSel(Atom(11, std::nullopt, mol2())).bonded().indices(), ElementsAre(6, 11));
}

TEST_F(Mol2MolfileReaderTest, GuessedBond)
{
    size_t const num_atoms = mol2()->properties().size<Atom>();
    for (size_t i = 0; i < num_atoms; i++)
    {
        for (auto bond : Atom(i, std::nullopt, mol2()).bonds())
        {
            EXPECT_FALSE(bond->guessed());
        }
    }
}

TEST_F(Mol2MolfileReaderTest, BondOrder)
{
    EXPECT_EQ(Atom(0, std::nullopt, mol2()).bond(1)->order(), 1);
    EXPECT_EQ(Atom(1, std::nullopt, mol2()).bond(2)->order(), 0);
    EXPECT_EQ(Atom(1, std::nullopt, mol2()).bond(3)->order(), 0);
    EXPECT_EQ(Atom(2, std::nullopt, mol2()).bond(4)->order(), 0);
    EXPECT_EQ(Atom(2, std::nullopt, mol2()).bond(7)->order(), 1);
    EXPECT_EQ(Atom(3, std::nullopt, mol2()).bond(5)->order(), 0);
    EXPECT_EQ(Atom(3, std::nullopt, mol2()).bond(8)->order(), 1);
    EXPECT_EQ(Atom(4, std::nullopt, mol2()).bond(6)->order(), 0);
    EXPECT_EQ(Atom(4, std::nullopt, mol2()).bond(9)->order(), 1);
    EXPECT_EQ(Atom(5, std::nullopt, mol2()).bond(6)->order(), 0);
    EXPECT_EQ(Atom(5, std::nullopt, mol2()).bond(10)->order(), 1);
    EXPECT_EQ(Atom(6, std::nullopt, mol2()).bond(11)->order(), 1);
}

TEST_F(Mol2MolfileReaderTest, GuessedOrder)
{
    size_t const num_atoms = mol2()->properties().size<Atom>();
    for (size_t i = 0; i < num_atoms; i++)
    {
        for (auto bond : Atom(i, std::nullopt, mol2()).bonds())
        {
            EXPECT_FALSE(bond->guessed_order()) << bond->atom1() << "-" << bond->atom2();
        }
    }
}

class PSFMolfileReaderTest : public BaseMolfileReaderTest
{
protected:
    void SetUp() override
    {
        ASSERT_THAT(psf(), NotNull());
        ASSERT_EQ(psf()->properties().size<Atom>(), 22);
    }

    MolData* psf()
    {
        static std::unique_ptr<MolData> data = read("dipeptide.psf", ".psf");
        return data.get();
    }
};

TEST_F(PSFMolfileReaderTest, AtomName)
{
    test_property<Atom, Name, 22>({"N", "CA", "C", "O", "N", "CD", "CA", "CB", "CG", "C", "O", "N", "CA", "C", "O", "N", "CD", "CA", "CB", "CG", "C", "O"}, psf());
}

TEST_F(PSFMolfileReaderTest, AtomType)
{
    test_property<Atom, Type, 22>({"NH1", "CT2", "C", "O", "N", "CP3", "CP1", "CP2", "CP2", "C", "O", "NH1", "CT2", "C", "O", "N", "CP3", "CP1", "CP2", "CP2", "C", "O"}, psf());
}

TEST_F(PSFMolfileReaderTest, Occupancy)
{
    test_no_property<Atom, Occupancy>(psf());
}

TEST_F(PSFMolfileReaderTest, TemperatureFactor)
{
    test_no_property<Atom, TemperatureFactor>(psf());
}

TEST_F(PSFMolfileReaderTest, Mass)
{
    test_property<Atom, Mass, 22>({14.007, 12.011, 12.011, 15.999, 14.007, 12.011, 12.011, 12.011, 12.011, 12.011, 15.999, 14.007, 12.011, 12.011, 15.999, 14.007, 12.011, 12.011, 12.011, 12.011, 12.011, 15.999}, psf());
}

TEST_F(PSFMolfileReaderTest, Radius)
{
    test_no_property<Atom, Radius>(psf());
}

TEST_F(PSFMolfileReaderTest, AtomicNumber)
{
    test_no_property<Atom, AtomicNumber>(psf());
}

TEST_F(PSFMolfileReaderTest, Charge)
{
    test_property<Atom, Charge, 22>({-0.47, -0.02, 0.51, -0.51, -0.29, 0.0, 0.02, -0.18, -0.18, 0.51, -0.51, -0.47, -0.02, 0.51, -0.51, -0.29, 0.0, 0.02, -0.18, -0.18, 0.51, -0.51}, psf());
}

TEST_F(PSFMolfileReaderTest, AlternateLocation)
{
    test_no_property<Atom, AlternateLocation>(psf());
}

TEST_F(PSFMolfileReaderTest, InsertionCode)
{
    test_property<Atom, InsertionCode, 22>({" ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " "}, psf());
}

TEST_F(PSFMolfileReaderTest, ResID)
{
    test_property<Residue, ResID, 4>({129, 130, 129, 130}, psf());
}

TEST_F(PSFMolfileReaderTest, ResName)
{
    test_property<Residue, ResName, 4>({"GLY", "PRO", "GLY", "PRO"}, psf());
}

TEST_F(PSFMolfileReaderTest, Bonded)
{
    EXPECT_THAT(AtomSel(Atom(0, std::nullopt, psf())).bonded().indices(), ElementsAre(0, 1));
    EXPECT_THAT(AtomSel(Atom(1, std::nullopt, psf())).bonded().indices(), ElementsAre(0, 1, 2));
    EXPECT_THAT(AtomSel(Atom(2, std::nullopt, psf())).bonded().indices(), ElementsAre(1, 2, 3, 4));
    EXPECT_THAT(AtomSel(Atom(3, std::nullopt, psf())).bonded().indices(), ElementsAre(2, 3));
    EXPECT_THAT(AtomSel(Atom(4, std::nullopt, psf())).bonded().indices(), ElementsAre(2, 4, 5, 6));
    EXPECT_THAT(AtomSel(Atom(5, std::nullopt, psf())).bonded().indices(), ElementsAre(4, 5, 8));
    EXPECT_THAT(AtomSel(Atom(6, std::nullopt, psf())).bonded().indices(), ElementsAre(4, 6, 7, 9));
    EXPECT_THAT(AtomSel(Atom(7, std::nullopt, psf())).bonded().indices(), ElementsAre(6, 7, 8));
    EXPECT_THAT(AtomSel(Atom(8, std::nullopt, psf())).bonded().indices(), ElementsAre(5, 7, 8));
    EXPECT_THAT(AtomSel(Atom(9, std::nullopt, psf())).bonded().indices(), ElementsAre(6, 9, 10));
    EXPECT_THAT(AtomSel(Atom(10, std::nullopt, psf())).bonded().indices(), ElementsAre(9, 10));
    EXPECT_THAT(AtomSel(Atom(11, std::nullopt, psf())).bonded().indices(), ElementsAre(11, 12));
    EXPECT_THAT(AtomSel(Atom(12, std::nullopt, psf())).bonded().indices(), ElementsAre(11, 12, 13));
    EXPECT_THAT(AtomSel(Atom(13, std::nullopt, psf())).bonded().indices(), ElementsAre(12, 13, 14, 15));
    EXPECT_THAT(AtomSel(Atom(14, std::nullopt, psf())).bonded().indices(), ElementsAre(13, 14));
    EXPECT_THAT(AtomSel(Atom(15, std::nullopt, psf())).bonded().indices(), ElementsAre(13, 15, 16, 17));
    EXPECT_THAT(AtomSel(Atom(16, std::nullopt, psf())).bonded().indices(), ElementsAre(15, 16, 19));
    EXPECT_THAT(AtomSel(Atom(17, std::nullopt, psf())).bonded().indices(), ElementsAre(15, 17, 18, 20));
    EXPECT_THAT(AtomSel(Atom(18, std::nullopt, psf())).bonded().indices(), ElementsAre(17, 18, 19));
    EXPECT_THAT(AtomSel(Atom(19, std::nullopt, psf())).bonded().indices(), ElementsAre(16, 18, 19));
    EXPECT_THAT(AtomSel(Atom(20, std::nullopt, psf())).bonded().indices(), ElementsAre(17, 20, 21));
    EXPECT_THAT(AtomSel(Atom(21, std::nullopt, psf())).bonded().indices(), ElementsAre(20, 21));
}

TEST_F(PSFMolfileReaderTest, GuessedBond)
{
    size_t const num_atoms = psf()->properties().size<Atom>();
    for (size_t i = 0; i < num_atoms; i++)
    {
        for (auto bond : Atom(i, std::nullopt, psf()).bonds())
        {
            EXPECT_FALSE(bond->guessed());
        }
    }
}

TEST_F(PSFMolfileReaderTest, BondOrder)
{
    size_t const num_atoms = psf()->properties().size<Atom>();
    for (size_t i = 0; i < num_atoms; i++)
    {
        for (auto bond : Atom(i, std::nullopt, psf()).bonds())
        {
            EXPECT_EQ(bond->order(), 1) << bond->atom1() << "-" << bond->atom2();
        }
    }
}

TEST_F(PSFMolfileReaderTest, GuessedOrder)
{
    size_t const num_atoms = psf()->properties().size<Atom>();
    for (size_t i = 0; i < num_atoms; i++)
    {
        for (auto bond : Atom(i, std::nullopt, psf()).bonds())
        {
            EXPECT_TRUE(bond->guessed_order()) << bond->atom1() << "-" << bond->atom2();
        }
    }
}

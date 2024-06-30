#include "utils.hpp"

#include <molpp/Atom.hpp>
#include <molpp/Residue.hpp>
#include <molpp/Chain.hpp>
#include <molpp/Segment.hpp>
#include <molpp/internal/MolData.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <array>

using namespace mol;
using namespace mol::internal;
using namespace testing;

//! Test fixture for create_moldata
class CreateMolDataTest : public ::testing::Test
{
public:
    CreateMolDataTest()
    : data(create_moldata(3, 2, 2, 2, 2))
    {
    }

    mol::internal::MolData data;
};

// create_moldata should create MolData with correct sizes
TEST_F(CreateMolDataTest, sizes)
{
    ASSERT_EQ(data.size<Atom>(), 6);
    ASSERT_EQ(data.size<Residue>(), 3);
    ASSERT_EQ(data.size<Chain>(), 2);
    ASSERT_EQ(data.size<Segment>(), 2);
}

// create_moldata should create MolData with correct number of frames
TEST_F(CreateMolDataTest, number_of_frames)
{
    size_t const num_frames = data.trajectory().num_frames();

    EXPECT_EQ(num_frames, 2);
}

// create_moldata should set correct values for atom properties
TEST_F(CreateMolDataTest, atom_properties)
{
    test_property(&Atom::atomic_number, std::to_array({0, 1, 2, 3, 4, 5}), data);
    test_property(&Atom::occupancy, std::to_array({0, 1, 2, 3, 4, 5}), data);
    test_property(&Atom::temperature_factor, std::to_array({0, 1, 2, 3, 4, 5}), data);
    test_property(&Atom::mass, std::to_array({0, 1, 2, 3, 4, 5}), data);
    test_property(&Atom::charge, std::to_array({0, 1, 2, 3, 4, 5}), data);
    test_property(&Atom::radius, std::to_array({0, 1, 2, 3, 4, 5}), data);
    test_property(&Atom::name, std::to_array({"A", "B", "C", "D", "E", "F"}), data);
    test_property(&Atom::type, std::to_array({"A", "B", "C", "D", "E", "F"}), data);
    test_property(&Atom::alternate_location, std::to_array({"A", "B", "C", "D", "E", "F"}), data);
    test_property(&Atom::insertion_code, std::to_array({"A", "B", "C", "D", "E", "F"}), data);
}

// create_moldata should set correct values for residue properties
TEST_F(CreateMolDataTest, residue_properties)
{
    test_property(&Residue::id, std::to_array({0, 1, 2}), data);
    test_property(&Residue::name, std::to_array({"A", "B", "C"}), data);
}

// create_moldata should set correct values for chain properties
TEST_F(CreateMolDataTest, chain_properties)
{
    test_property(&Chain::name, std::to_array({"A", "B"}), data);
}

// create_moldata should set correct values for segment properties
TEST_F(CreateMolDataTest, segment_properties)
{
    test_property(&Segment::name, std::to_array({"A", "B"}), data);
}

// create_moldata should link atoms to residues
TEST_F(CreateMolDataTest, link_atoms_to_residues)
{
    test_property(&Atom::residue_index, std::to_array({0, 0, 1, 1, 2, 2}), data);
}

// create_moldata should link residues to chains
TEST_F(CreateMolDataTest, link_residues_to_chains)
{
    test_property(&Residue::chain_index, std::to_array({0, 1, 0}), data);
}

// create_moldata should link residues to segments
TEST_F(CreateMolDataTest, link_residues_to_segments)
{
    test_property(&Residue::segment_index, std::to_array({0, 1, 0}), data);
}

// create_moldata should add bonds only between first atoms in each residue
TEST_F(CreateMolDataTest, bonds)
{
    BondData const& bond_data = data.bonds();

    EXPECT_THAT(bond_data.bonded(0), UnorderedElementsAre(0, 2));
    EXPECT_THAT(bond_data.bonded(1), UnorderedElementsAre());
    EXPECT_THAT(bond_data.bonded(2), UnorderedElementsAre(0, 2, 4));
    EXPECT_THAT(bond_data.bonded(3), UnorderedElementsAre());
    EXPECT_THAT(bond_data.bonded(4), UnorderedElementsAre(2, 4));
    EXPECT_THAT(bond_data.bonded(5), UnorderedElementsAre());
}

// create_moldata should set correct coordinates for each frame
TEST_F(CreateMolDataTest, positions)
{
    Trajectory const& trajectory = data.trajectory();

    ASSERT_EQ(trajectory.num_frames(), 2);
    EXPECT_THAT(trajectory.timestep(0).coords().reshaped(), ElementsAre(0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5));
    EXPECT_THAT(trajectory.timestep(1).coords().reshaped(), ElementsAre(0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5));
}
